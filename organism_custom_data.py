# ------------------
# Script to retrieve PDB structures for a custom gene list
# ------------------
import pandas as pd
import os
import requests
from typing import List, Optional, Dict, Tuple
import time
import re
import logging

# Configuration
MAX_RETRIES = 3
REQUEST_DELAY = 0.5

def setup_logging(organism_path: str) -> logging.Logger:
    '''Set up logging for the script'''
    log_file = os.path.join(organism_path, "pdb_retrieval.log")
    
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    
    # Clear any existing handlers
    logger.handlers.clear()
    
    # File handler
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setLevel(logging.INFO)
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    
    # Formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)
    
    # Add handlers
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return logger

def parse_pdb_data_from_uniprot(uniprot_text: str) -> List[Tuple[str, float, str]]:
    '''Parse PDB IDs, resolutions, and methods from UniProt text data'''
    pdb_data = []
    lines = uniprot_text.split('\n')
    
    for line in lines:
        if line.startswith("DR   PDB;") or "PDB;" in line and line.startswith("DR"):
            parts = line.strip().split(';')
            if len(parts) >= 4:
                pdb_id = parts[1].strip()
                method = parts[2].strip()
                resolution_str = parts[3].strip()
                
                resolution_match = re.search(r'([\d\.]+)\s*A', resolution_str)
                if resolution_match:
                    try:
                        resolution = float(resolution_match.group(1))
                        pdb_data.append((pdb_id, resolution, method))
                    except ValueError:
                        continue
                elif resolution_str == "-":
                    pdb_data.append((pdb_id, float('inf'), method))
                elif resolution_str in ["NMR", "Model", "THEORETICAL MODEL"]:
                    pdb_data.append((pdb_id, float('inf'), method))
    
    return pdb_data

def get_best_pdb_structure(pdb_data: List[Tuple[str, float, str]]) -> Tuple[Optional[str], Optional[float], Optional[str]]:
    '''Get the PDB ID with the best (lowest) resolution'''
    if not pdb_data:
        return None, None, None
    
    sorted_data = sorted(pdb_data, key=lambda x: x[1])
    best_pdb_id, best_resolution, best_method = sorted_data[0]
    return best_pdb_id, best_resolution, best_method

def download_pdb_file(pdb_id: str, output_dir: str, logger: logging.Logger) -> Optional[str]:
    '''Download PDB file from RCSB PDB'''
    # Try multiple URL formats
    url_formats = [
        f"https://files.rcsb.org/download/{pdb_id.lower()}.pdb",
        f"https://files.rcsb.org/view/{pdb_id}.pdb",
    ]
    
    for url in url_formats:
        try:
            response = requests.get(url, timeout=30)
            
            if response.status_code == 200:
                content = response.text
                if content and ('HEADER' in content or 'ATOM' in content or 'REMARK' in content):
                    pdb_file_path = os.path.join(output_dir, f"{pdb_id}.pdb")
                    
                    with open(pdb_file_path, "w") as f:
                        f.write(content)
                    
                    logger.info(f"  Downloaded PDB file: {pdb_file_path} from {url}")
                    print(f"  Downloaded PDB file: {pdb_file_path}")
                    return pdb_file_path
                else:
                    logger.warning(f"  Response for {pdb_id} from {url} doesn't look like a valid PDB file")
            elif response.status_code == 404:
                logger.debug(f"  PDB {pdb_id} not found at {url}")
                continue
            else:
                logger.debug(f"  HTTP {response.status_code} for {pdb_id} from {url}")
                
        except requests.exceptions.RequestException as e:
            logger.debug(f"  Error downloading PDB {pdb_id} from {url}: {e}")
            continue
    
    logger.warning(f"  Failed to download PDB {pdb_id} from all available sources")
    print(f"  Failed to download PDB {pdb_id}")
    return None

def download_pdb_file_with_retry(pdb_id: str, output_dir: str, logger: logging.Logger, max_retries: int = 3) -> Optional[str]:
    '''Download PDB file with retry logic'''
    for attempt in range(max_retries):
        # CORRECTED: Only pass 3 arguments to download_pdb_file
        result = download_pdb_file(pdb_id, output_dir, logger)
        if result:
            if attempt > 0:
                logger.info(f"    Successfully downloaded {pdb_id} on attempt {attempt+1}")
            return result
        
        if attempt < max_retries - 1:
            wait_time = 2 ** attempt
            logger.debug(f"    Waiting {wait_time} seconds before retry...")
            time.sleep(wait_time)
    
    logger.warning(f"    Failed to download PDB {pdb_id} after {max_retries} attempts")
    return None

def search_uniprot_for_gene(gene_name: str, organism: str, logger: logging.Logger) -> Optional[str]:
    '''
    Search UniProt for a gene name and return the first UniProt accession
    
    :param gene_name: Gene name to search
    :type gene_name: str
    :param organism: Organism name to filter results
    :type organism: str
    :param logger: Logger object
    :type logger: logging.Logger
    :return: UniProt accession or None
    :rtype: str or None
    '''
    try:
        # UniProt search API
        search_url = "https://rest.uniprot.org/uniprotkb/search"
        
        # Search for the gene name with organism filter
        query = f"gene:{gene_name} AND organism_name:\"{organism}\""
        
        params = {
            "query": query,
            "format": "tsv",
            "fields": "accession,gene_names,protein_name,organism_name",
            "size": 5  # Get top 5 results
        }
        
        logger.info(f"Searching UniProt for gene: {gene_name} in organism: {organism}")
        response = requests.get(search_url, params=params, timeout=30)
        response.raise_for_status()
        
        # Parse the TSV response
        lines = response.text.strip().split('\n')
        if len(lines) <= 1:  # Only header or empty
            logger.warning(f"No results found for gene: {gene_name} in organism: {organism}")
            
            # Try a broader search without organism filter
            logger.info(f"Trying broader search for gene: {gene_name} (any organism)")
            params["query"] = f"gene:{gene_name}"
            response = requests.get(search_url, params=params, timeout=30)
            response.raise_for_status()
            
            lines = response.text.strip().split('\n')
            if len(lines) <= 1:
                logger.warning(f"No results found for gene: {gene_name} in any organism")
                return None
        
        # Get the first result (skip header)
        first_result = lines[1].split('\t')
        if len(first_result) > 0:
            uniprot_accession = first_result[0]
            found_organism = first_result[3] if len(first_result) > 3 else "Unknown"
            logger.info(f"Found UniProt accession {uniprot_accession} for gene {gene_name} (Organism: {found_organism})")
            return uniprot_accession
        
        return None
        
    except requests.exceptions.RequestException as e:
        logger.error(f"Error searching UniProt for {gene_name}: {e}")
        return None
    except Exception as e:
        logger.error(f"Unexpected error searching for {gene_name}: {e}")
        return None

def retrieve_uniprot_data(uniprot_id: str, output_dir: str, logger: logging.Logger) -> Optional[str]:
    '''
    Retrieve UniProt data for a given accession
    
    :param uniprot_id: UniProt accession
    :type uniprot_id: str
    :param output_dir: Directory to save UniProt data
    :type output_dir: str
    :param logger: Logger object
    :type logger: logging.Logger
    :return: Path to UniProt file or None
    :rtype: str or None
    '''
    try:
        uniprot_api = "https://rest.uniprot.org/uniprotkb/"
        uniprot_url = f"{uniprot_api}{uniprot_id}"
        
        logger.info(f"Retrieving UniProt data for: {uniprot_id}")
        response = requests.get(uniprot_url, params={"format": "txt"}, timeout=30)
        response.raise_for_status()
        
        # Save UniProt data
        uniprot_file = os.path.join(output_dir, f"{uniprot_id}_uniprot.txt")
        with open(uniprot_file, "w") as f:
            f.write(response.text)
        
        logger.info(f"Saved UniProt data to: {uniprot_file}")
        return uniprot_file
        
    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 404:
            logger.error(f"UniProt accession not found: {uniprot_id}")
        else:
            logger.error(f"HTTP error retrieving UniProt data for {uniprot_id}: {e}")
        return None
    except requests.exceptions.RequestException as e:
        logger.error(f"Network error retrieving UniProt data for {uniprot_id}: {e}")
        return None
    except Exception as e:
        logger.error(f"Unexpected error retrieving UniProt data for {uniprot_id}: {e}")
        return None

def process_gene_list(
    gene_list_path: str,
    organism: str,
    output_path: str = "custom_genes/"
) -> Optional[pd.DataFrame]:
    '''
    Process a custom gene list and retrieve PDB structures
    
    :param gene_list_path: Path to text file with gene names (one per line)
    :type gene_list_path: str
    :param organism: Organism name for UniProt search
    :type organism: str
    :param output_path: Directory to store results
    :type output_path: str
    :return: DataFrame with results
    :rtype: pd.DataFrame or None
    '''
    
    # 1. Verify if input file exists
    if not os.path.exists(gene_list_path):
        print(f"Gene list file not found at: {gene_list_path}")
        return None
    
    # 2. Read gene list
    try:
        with open(gene_list_path, 'r') as f:
            gene_list = [line.strip() for line in f if line.strip()]
        
        if not gene_list:
            print("Gene list file is empty")
            return None
        
        print(f"Loaded {len(gene_list)} genes from: {gene_list_path}")
        print(f"Target organism: {organism}")
        
    except Exception as e:
        print(f"Error reading gene list file: {e}")
        return None
    
    # 3. Create output directory
    if not os.path.exists(output_path):
        try:
            os.makedirs(output_path)
            print(f"Created output directory: {output_path}")
        except Exception as e:
            print(f"Error creating output directory: {e}")
            return None
    
    # 4. Setup logging
    logger = setup_logging(output_path)
    logger.info(f"Starting processing of custom gene list: {gene_list_path}")
    logger.info(f"Target organism: {organism}")
    logger.info(f"Number of genes to process: {len(gene_list)}")
    
    # 5. Process each gene
    results = []
    summary_data = []
    
    for idx, gene_name in enumerate(gene_list):
        logger.info(f"\nProcessing {idx+1}/{len(gene_list)}: {gene_name}")
        print(f"\nProcessing {idx+1}/{len(gene_list)}: {gene_name}")
        
        # Add delay to avoid overloading the API
        if idx > 0:
            time.sleep(REQUEST_DELAY)
        
        # Search for UniProt accession with organism filter
        uniprot_id = search_uniprot_for_gene(gene_name, organism, logger)
        
        if not uniprot_id:
            results.append({
                'gene': gene_name,
                'organism': organism,
                'uniprot': 'Not Found',
                'pdb': 'Not Found',
                'resolution': 'Not Found',
                'method': 'Not Found',
                'pdb_file': 'Not Found',
                'uniprot_file': 'Not Found',
                'all_pdbs': 'Not Found'
            })
            continue
        
        # Retrieve UniProt data
        uniprot_file = retrieve_uniprot_data(uniprot_id, output_path, logger)
        
        if not uniprot_file:
            results.append({
                'gene': gene_name,
                'organism': organism,
                'uniprot': uniprot_id,
                'pdb': 'UniProt Error',
                'resolution': 'UniProt Error',
                'method': 'UniProt Error',
                'pdb_file': 'UniProt Error',
                'uniprot_file': 'UniProt Error',
                'all_pdbs': 'UniProt Error'
            })
            continue
        
        # Read UniProt data to parse PDB information
        try:
            with open(uniprot_file, 'r') as f:
                uniprot_text = f.read()
            
            # Parse PDB IDs and resolutions
            pdb_data = parse_pdb_data_from_uniprot(uniprot_text)
            
            if not pdb_data:
                logger.info(f"No PDB structures found for {uniprot_id}")
                print(f"  No PDB structures found")
                
                results.append({
                    'gene': gene_name,
                    'organism': organism,
                    'uniprot': uniprot_id,
                    'pdb': 'None',
                    'resolution': 'None',
                    'method': 'None',
                    'pdb_file': 'None',
                    'uniprot_file': uniprot_file,
                    'all_pdbs': 'None'
                })
                continue
            
            logger.info(f"Found {len(pdb_data)} PDB structures for {uniprot_id}")
            print(f"  Found {len(pdb_data)} PDB structures")
            
            # Get best PDB structure
            best_pdb_id, best_resolution, best_method = get_best_pdb_structure(pdb_data)
            
            # Download PDB file
            pdb_file_path = None
            downloaded_pdb_id = None
            downloaded_resolution = None
            downloaded_method = None
            
            if best_pdb_id:
                # Sort PDB data by resolution
                sorted_pdb_data = sorted(pdb_data, key=lambda x: x[1])
                
                # Try each PDB in order of best resolution
                for pdb_id, resolution, method in sorted_pdb_data:
                    logger.info(f"  Attempting to download {pdb_id} ({resolution} Å, {method})")
                    
                    pdb_file_path = download_pdb_file_with_retry(pdb_id, output_path, logger, MAX_RETRIES)
                    if pdb_file_path:
                        downloaded_pdb_id = pdb_id
                        downloaded_resolution = resolution
                        downloaded_method = method
                        break
            
            # Store results
            result_entry = {
                'gene': gene_name,
                'organism': organism,
                'uniprot': uniprot_id,
                'pdb': downloaded_pdb_id if downloaded_pdb_id else 'None',
                'resolution': downloaded_resolution if downloaded_pdb_id else 'None',
                'method': downloaded_method if downloaded_pdb_id else 'None',
                'pdb_file': pdb_file_path,
                'uniprot_file': uniprot_file,
                'all_pdbs': ', '.join([f"{data[0]}({data[1]}Å)" for data in pdb_data]) if pdb_data else 'None'
            }
            
            results.append(result_entry)
            
            # Store summary data if PDB was downloaded
            if downloaded_pdb_id:
                summary_data.append({
                    'gene': gene_name,
                    'organism': organism,
                    'uniprot': uniprot_id,
                    'pdb': downloaded_pdb_id,
                    'resolution': downloaded_resolution,
                    'method': downloaded_method
                })
            
            logger.info(f"  Result: PDB={downloaded_pdb_id}, Resolution={downloaded_resolution}Å, Method={downloaded_method}")
            
        except Exception as e:
            error_msg = f"Error processing UniProt data for {uniprot_id}: {e}"
            logger.error(error_msg)
            print(f"Error: {error_msg}")
            
            results.append({
                'gene': gene_name,
                'organism': organism,
                'uniprot': uniprot_id,
                'pdb': 'Processing Error',
                'resolution': 'Processing Error',
                'method': 'Processing Error',
                'pdb_file': 'Processing Error',
                'uniprot_file': uniprot_file if 'uniprot_file' in locals() else 'Error',
                'all_pdbs': 'Processing Error'
            })
    
    # 6. Create results DataFrame
    results_df = pd.DataFrame(results)
    
    # 7. Save detailed results with dynamic naming
    # Get the base name of the gene list file without extension
    gene_list_basename = os.path.basename(gene_list_path)
    gene_list_name = os.path.splitext(gene_list_basename)[0]
    
    detailed_csv = os.path.join(output_path, f"{gene_list_name}_results.csv")
    results_df.to_csv(detailed_csv, index=False)
    print(f"\nDetailed results saved to: {detailed_csv}")
    logger.info(f"Detailed results saved to: {detailed_csv}")
    
    # 8. Save summary results with dynamic naming
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_csv = os.path.join(output_path, f"{gene_list_name}_summary.csv")
        summary_df.to_csv(summary_csv, index=False)
        print(f"Summary results saved to: {summary_csv}")
        logger.info(f"Summary results saved to: {summary_csv}")
    else:
        summary_df = pd.DataFrame()
        logger.warning("No PDB files were successfully downloaded")
        print("Warning: No PDB files were successfully downloaded")
    
    # 9. Print summary
    print(f"\n=== Summary ===")
    print(f"Organism: {organism}")
    print(f"Total genes processed: {len(results)}")
    
    # Count successful PDB retrievals
    successful_pdbs = results_df[
        (results_df['pdb'].notna()) & 
        (results_df['pdb'] != 'Error') & 
        (results_df['pdb'] != 'None') &
        (results_df['pdb'] != 'Not Found') &
        (results_df['pdb'] != 'UniProt Error') &
        (results_df['pdb'] != 'Processing Error')
    ].shape[0]
    
    # Count genes with UniProt found
    uniprot_found = results_df[
        (results_df['uniprot'].notna()) & 
        (results_df['uniprot'] != 'Not Found')
    ].shape[0]
    
    print(f"Genes with UniProt accession found: {uniprot_found}")
    print(f"PDB files successfully downloaded: {successful_pdbs}")
    print(f"Output directory: {output_path}")
    
    logger.info(f"=== Summary ===")
    logger.info(f"Organism: {organism}")
    logger.info(f"Total genes processed: {len(results)}")
    logger.info(f"Genes with UniProt accession found: {uniprot_found}")
    logger.info(f"PDB files successfully downloaded: {successful_pdbs}")
    
    # Calculate resolution statistics
    if successful_pdbs > 0:
        valid_res_df = results_df[
            (results_df['resolution'].notna()) & 
            (results_df['resolution'] != 'Error') & 
            (results_df['resolution'] != 'None') &
            (results_df['resolution'] != 'Not Found') &
            (results_df['resolution'] != 'UniProt Error') &
            (results_df['resolution'] != 'Processing Error')
        ].copy()
        
        valid_res_df['resolution_numeric'] = pd.to_numeric(valid_res_df['resolution'], errors='coerce')
        valid_res_df = valid_res_df.dropna(subset=['resolution_numeric'])
        
        if not valid_res_df.empty:
            avg_resolution = valid_res_df['resolution_numeric'].mean()
            min_resolution = valid_res_df['resolution_numeric'].min()
            max_resolution = valid_res_df['resolution_numeric'].max()
            
            print(f"Average resolution: {avg_resolution:.2f} Å")
            print(f"Best resolution: {min_resolution} Å")
            print(f"Worst resolution: {max_resolution} Å")
            
            logger.info(f"Average resolution: {avg_resolution:.2f} Å")
            logger.info(f"Best resolution: {min_resolution} Å")
            logger.info(f"Worst resolution: {max_resolution} Å")
    
    return results_df, summary_df if summary_data else None

def main():
    '''
    Main function to run the custom gene list script
    '''
    print("=== Custom Gene List PDB Retrieval Script ===\n")
    
    # Configuration
    gene_list_file = "genes/custom_genes.txt"  # Change this to your gene list file
    organism = "Escherichia coli K12"  # Change this to your target organism
    
    # Create output directory name that includes organism
    organism_safe = organism.lower().replace(" ", "_").replace("/", "_")
    output_directory = f"genes/{organism_safe}"
    
    print(f"Gene list file: {gene_list_file}")
    print(f"Target organism: {organism}")
    print(f"Output directory: {output_directory}")
    print(f"Maximum retries per PDB: {MAX_RETRIES}")
    print(f"Request delay: {REQUEST_DELAY} seconds\n")
    
    # Check if output directory already exists
    if os.path.exists(output_directory):
        print(f"Output directory already exists: {output_directory}")
        response = input("Do you want to overwrite? (yes/no): ").strip().lower()
        if response != 'yes':
            print("Exiting script.")
            return
    
    # Process the gene list
    results, summary = process_gene_list(
        gene_list_path=gene_list_file,
        organism=organism,
        output_path=output_directory
    )
    
    if results is not None:
        print(f"\n{'='*60}")
        print("Script completed successfully!")
        print(f"\nResults saved in: {output_directory}")
        
        # Show outcome distribution
        print(f"\nOutcome distribution:")
        outcome_counts = results['pdb'].value_counts()
        for outcome, count in outcome_counts.items():
            percentage = (count / len(results)) * 100
            print(f"  {outcome}: {count} ({percentage:.1f}%)")
    
    else:
        print(f"\nScript failed to complete.")

if __name__ == "__main__":
    main()