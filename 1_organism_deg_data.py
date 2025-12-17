# ------------------
# Script to retrieve gene list from organism
# ------------------
import pandas as pd
import os
import requests
from typing import List, Optional, Dict, Tuple
import time
import re
import logging
from datetime import datetime

# FUNCTIONS ----------------------------
def setup_logging(organism_path: str) -> logging.Logger:
    '''
    Set up logging for the script
    
    :param organism_path: Path to organism directory
    :type organism_path: str
    :return: Logger object
    :rtype: logging.Logger
    '''
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

def retrieve_organism_data(
    organism: str,
    deg_path: str = "genes/deg.csv",
    output_path: str = "genes/"
) -> Optional[pd.DataFrame]:
    '''
    Retrieve gene data for a specific organism from DEG database and get PDB files
    
    :param organism: Name of organism from DEG Database
    :type organism: str
    :param deg_path: Path to DEG database CSV file
    :type deg_path: str
    :param output_path: Directory to store PDB files and organism data
    :type output_path: str
    :return: DataFrame with organism genes and PDB information
    :rtype: pd.DataFrame or None
    '''
    
    # 1. Verify if output directory exists, create if not
    if not os.path.exists(output_path):
        try:
            os.makedirs(output_path)
            print(f"Created directory: {output_path}")
        except Exception as e:
            print(f"Error creating directory {output_path}: {e}")
            return None
    
    # 2. Verify if Database exists
    if not os.path.exists(deg_path):
        print(f"Database not found at: {deg_path}")
        return None
    
    # 3. Open Database
    try:
        deg_df = pd.read_csv(deg_path)
        print(f"Successfully loaded database with {len(deg_df)} entries")
    except Exception as e:
        print(f"Error loading database: {e}")
        return None
    
    # 4. Filter Database by organism
    try:
        # Correct the filtering syntax
        organism_df = deg_df[deg_df["organism"] == organism].copy()
        
        if organism_df.empty:
            print(f"Organism '{organism}' not found in database.")
            print(f"Available organisms: {deg_df['organism'].unique()[:10]}...")
            return None
        
        print(f"Found {len(organism_df)} entries for organism: {organism}")
        
        # Get UniProt IDs and gene names
        if "uniprot" not in organism_df.columns:
            print("'uniprot' column not found in database.")
            return None
        
        # Check if 'gene' column exists, if not use 'locus' or create placeholder
        gene_column = None
        for col in ['gene', 'locus', 'gene_name']:
            if col in organism_df.columns:
                gene_column = col
                break
        
        if gene_column is None:
            print("No gene name column found. Creating placeholder.")
            organism_df['gene'] = organism_df['uniprot']
            gene_column = 'gene'
        
    except KeyError as e:
        print(f"Column error: {e}. Available columns: {deg_df.columns.tolist()}")
        return None
    except Exception as e:
        print(f"Error filtering data: {e}")
        return None
    
    # 5. Retrieve UniProt data and PDB information
    organism_path = os.path.join(output_path, organism.lower().replace(" ", "_"))
    if not os.path.exists(organism_path):
        os.makedirs(organism_path)
        print(f"Created organism directory: {organism_path}")
    
    # Setup logging
    logger = setup_logging(organism_path)
    logger.info(f"Starting PDB retrieval for organism: {organism}")
    logger.info(f"Total entries to process: {len(organism_df)}")
    
    # Lists to store results
    results = []
    all_uniprot_data = []
    
    # Process each entry
    for idx, row in organism_df.iterrows():
        uniprot_id = row["uniprot"]
        gene_name = row[gene_column] if pd.notna(row[gene_column]) else uniprot_id
        
        logger.info(f"Processing {idx+1}/{len(organism_df)}: {uniprot_id} (Gene: {gene_name})")
        print(f"Processing {idx+1}/{len(organism_df)}: {uniprot_id}")
        
        try:
            # Add delay to avoid overloading the API
            if idx > 0:
                time.sleep(0.5)
            
            # Retrieve UniProt data
            uniprot_api = "https://rest.uniprot.org/uniprotkb/"
            uniprot_url = f"{uniprot_api}{uniprot_id}"
            
            response = requests.get(uniprot_url, params={"format": "txt"})
            response.raise_for_status()
            
            # Save UniProt data
            uniprot_file = os.path.join(organism_path, f"{uniprot_id}_uniprot.txt")
            with open(uniprot_file, "w") as f:
                f.write(response.text)
            
            # Parse PDB IDs and resolutions from UniProt data
            pdb_data = parse_pdb_data_from_uniprot(response.text)
            
            # Get best PDB ID (lowest resolution)
            best_pdb_id, best_resolution, best_method = get_best_pdb_structure(pdb_data)
            
            # Download PDB file - try multiple structures if needed
            pdb_file_path = None
            downloaded_pdb_id = None
            
            if best_pdb_id:
                # Sort PDB data by resolution
                sorted_pdb_data = sorted(pdb_data, key=lambda x: x[1])
                
                # Try each PDB in order of best resolution
                for pdb_id, resolution, method in sorted_pdb_data:
                    pdb_file_path = download_pdb_file(pdb_id, organism_path, logger)
                    if pdb_file_path:
                        downloaded_pdb_id = pdb_id
                        # Update best resolution and method to match downloaded structure
                        best_resolution = resolution
                        best_method = method
                        break
            
            # Store results
            result_entry = {
                'gene': gene_name,
                'uniprot': uniprot_id,
                'pdb': downloaded_pdb_id if downloaded_pdb_id else 'None',
                'resolution': best_resolution if downloaded_pdb_id else 'None',
                'method': best_method if downloaded_pdb_id else 'None',
                'pdb_file': pdb_file_path,
                'uniprot_file': uniprot_file,
                'all_pdbs': ', '.join([f"{data[0]}({data[1]}Å)" for data in pdb_data]) if pdb_data else 'None'
            }
            
            results.append(result_entry)
            
            # Also store data for final summary
            if downloaded_pdb_id:
                all_uniprot_data.append({
                    'gene': gene_name,
                    'uniprot': uniprot_id,
                    'pdb': downloaded_pdb_id,
                    'resolution': best_resolution,
                    'method': best_method
                })
            
            logger.info(f"  Result: PDB={downloaded_pdb_id}, Resolution={best_resolution}Å, Method={best_method}")
            
        except requests.exceptions.RequestException as e:
            error_msg = f"Error retrieving data for {uniprot_id}: {e}"
            logger.error(error_msg)
            print(f"Error: {error_msg}")
            
            results.append({
                'gene': gene_name,
                'uniprot': uniprot_id,
                'pdb': 'Error',
                'resolution': 'Error',
                'method': 'Error',
                'pdb_file': 'Error',
                'uniprot_file': 'Error',
                'all_pdbs': 'Error'
            })
        except Exception as e:
            error_msg = f"Unexpected error for {uniprot_id}: {e}"
            logger.error(error_msg)
            print(f"Error: {error_msg}")
            
            results.append({
                'gene': gene_name,
                'uniprot': uniprot_id,
                'pdb': 'Error',
                'resolution': 'Error',
                'method': 'Error',
                'pdb_file': 'Error',
                'uniprot_file': 'Error',
                'all_pdbs': 'Error'
            })
    
    # 6. Create detailed results DataFrame
    results_df = pd.DataFrame(results)
    
    # 7. Save detailed results to CSV in organism directory
    output_csv = os.path.join(organism_path, f"{organism.lower().replace(' ', '_')}_detailed_results.csv")
    results_df.to_csv(output_csv, index=False)
    print(f"\nDetailed results saved to: {output_csv}")
    logger.info(f"Detailed results saved to: {output_csv}")
    
    # 8. Save final summary CSV in genes directory (not organism subdirectory)
    if all_uniprot_data:
        summary_df = pd.DataFrame(all_uniprot_data)
        summary_filename = f"{organism.lower().replace(' ', '_')}_summary.csv"
        summary_path = os.path.join(output_path, summary_filename)
        summary_df.to_csv(summary_path, index=False)
        print(f"Summary CSV saved to: {summary_path}")
        logger.info(f"Summary CSV saved to: {summary_path}")
    
    # 9. Print summary
    print(f"\n=== Summary ===")
    print(f"Organism: {organism}")
    print(f"Total entries processed: {len(results)}")
    
    # Count successful PDB retrievals
    successful_pdbs = results_df[
        (results_df['pdb'].notna()) & 
        (results_df['pdb'] != 'Error') & 
        (results_df['pdb'] != 'None')
    ].shape[0]
    
    print(f"PDB files successfully downloaded: {successful_pdbs}")
    print(f"Output directory: {organism_path}")
    
    logger.info(f"=== Summary ===")
    logger.info(f"Organism: {organism}")
    logger.info(f"Total entries processed: {len(results)}")
    logger.info(f"PDB files successfully downloaded: {successful_pdbs}")
    
    # Calculate statistics
    if successful_pdbs > 0:
        valid_res_df = results_df[
            (results_df['resolution'].notna()) & 
            (results_df['resolution'] != 'Error') & 
            (results_df['resolution'] != 'None')
        ].copy()
        
        # Convert resolution to numeric for calculations
        valid_res_df['resolution_numeric'] = pd.to_numeric(valid_res_df['resolution'], errors='coerce')
        
        avg_resolution = valid_res_df['resolution_numeric'].mean()
        min_resolution = valid_res_df['resolution_numeric'].min()
        max_resolution = valid_res_df['resolution_numeric'].max()
        
        print(f"Average resolution: {avg_resolution:.2f} Å")
        print(f"Best resolution: {min_resolution} Å")
        print(f"Worst resolution (among downloaded): {max_resolution} Å")
        
        logger.info(f"Average resolution: {avg_resolution:.2f} Å")
        logger.info(f"Best resolution: {min_resolution} Å")
        logger.info(f"Worst resolution: {max_resolution} Å")
    
    return results_df, summary_df if all_uniprot_data else None

def parse_pdb_data_from_uniprot(uniprot_text: str) -> List[Tuple[str, float, str]]:
    '''
    Parse PDB IDs, resolutions, and methods from UniProt text data
    
    :param uniprot_text: Text response from UniProt API
    :type uniprot_text: str
    :return: List of tuples (PDB ID, resolution, method)
    :rtype: List[Tuple[str, float, str]]
    '''
    pdb_data = []
    lines = uniprot_text.split('\n')
    
    for line in lines:
        if line.startswith("DR   PDB;") or "PDB;" in line and line.startswith("DR"):
            # Example line: "DR   PDB; 3J3V; EM; 13.30 A; U=1-103."
            
            # Extract PDB ID, method, and resolution
            parts = line.strip().split(';')
            if len(parts) >= 4:
                pdb_id = parts[1].strip()
                method = parts[2].strip()  # EM, X-RAY, NMR, etc.
                resolution_str = parts[3].strip()
                
                # Extract numerical resolution value
                resolution_match = re.search(r'([\d\.]+)\s*A', resolution_str)
                if resolution_match:
                    try:
                        resolution = float(resolution_match.group(1))
                        pdb_data.append((pdb_id, resolution, method))
                    except ValueError:
                        # Skip if resolution can't be converted to float
                        continue
                # Some entries might have "-" or other markers for no resolution
                elif resolution_str == "-":
                    # Assign a very high resolution as placeholder for sorting
                    pdb_data.append((pdb_id, float('inf'), method))
                # Some might be "NMR" or other non-numerical entries
                elif resolution_str in ["NMR", "Model"]:
                    pdb_data.append((pdb_id, float('inf'), method))
    
    return pdb_data

def get_best_pdb_structure(pdb_data: List[Tuple[str, float, str]]) -> Tuple[Optional[str], Optional[float], Optional[str]]:
    '''
    Get the PDB ID with the best (lowest) resolution
    
    :param pdb_data: List of tuples (PDB ID, resolution, method)
    :type pdb_data: List[Tuple[str, float, str]]
    :return: Best PDB ID, resolution, and method
    :rtype: Tuple[Optional[str], Optional[float], Optional[str]]
    '''
    if not pdb_data:
        return None, None, None
    
    # Sort by resolution (lower is better)
    # Handle infinite resolutions (those marked with "-") by sorting them last
    sorted_data = sorted(pdb_data, key=lambda x: x[1])
    
    # Get the best (lowest resolution) structure
    best_pdb_id, best_resolution, best_method = sorted_data[0]
    
    return best_pdb_id, best_resolution, best_method

def download_pdb_file(pdb_id: str, output_dir: str, logger: logging.Logger) -> Optional[str]:
    '''
    Download PDB file from RCSB PDB
    
    :param pdb_id: PDB ID (e.g., '1ABC')
    :type pdb_id: str
    :param output_dir: Directory to save PDB file
    :type output_dir: str
    :param logger: Logger object for logging
    :type logger: logging.Logger
    :return: Path to downloaded file or None
    :rtype: str or None
    '''
    # Try multiple URL formats
    url_formats = [
        f"https://files.rcsb.org/download/{pdb_id.lower()}.pdb",
        f"https://files.rcsb.org/view/{pdb_id}.pdb",
        f"https://www.rcsb.org/structure/{pdb_id}",  # This might redirect to download
    ]
    
    for url in url_formats:
        time.sleep(0.5)
        try:
            response = requests.get(url, timeout=30)
            
            if response.status_code == 200:
                # Check if the response is actually a PDB file
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
            else:
                logger.debug(f"  Failed to download {pdb_id} from {url}: HTTP {response.status_code}")
                
        except requests.exceptions.RequestException as e:
            logger.debug(f"  Error downloading PDB {pdb_id} from {url}: {e}")
            continue
    
    logger.warning(f"  Failed to download PDB {pdb_id} from all available sources")
    print(f"  Failed to download PDB {pdb_id}")
    return None

def merge_all_summaries(output_path: str = "genes/"):
    '''
    Merge all organism summary files into one master CSV
    
    :param output_path: Directory containing summary files
    :type output_path: str
    :return: Combined DataFrame
    :rtype: pd.DataFrame
    '''
    summary_files = [f for f in os.listdir(output_path) if f.endswith('_summary.csv')]
    
    all_data = []
    
    for summary_file in summary_files:
        file_path = os.path.join(output_path, summary_file)
        try:
            df = pd.read_csv(file_path)
            df['organism'] = summary_file.replace('_summary.csv', '').replace('_', ' ')
            all_data.append(df)
            print(f"Loaded: {summary_file} ({len(df)} entries)")
        except Exception as e:
            print(f"Error loading {summary_file}: {e}")
    
    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)
        master_file = os.path.join(output_path, "all_organisms_summary.csv")
        combined_df.to_csv(master_file, index=False)
        print(f"\nMaster summary saved to: {master_file}")
        print(f"Total entries: {len(combined_df)}")
        return combined_df
    else:
        print("No summary files found")
        return pd.DataFrame()

def main():
    '''
    Main function to run the script for all organisms
    '''
    deg_path = "genes/deg.csv"  # Path to Database
    organisms = ["Bacillus subtilis 168"]

    print("=== Gene Retrieval Script ===\n")

    all_detailed_results = {}
    all_summaries = []

    for organism in organisms:
        print(f"\n{'='*50}")
        print(f"Processing organism: {organism}")
        print(f"{'='*50}")

        # Check if organism folder already exists
        organism_folder = os.path.join("genes", organism.lower().replace(" ", "_"))
        if os.path.exists(organism_folder):
            print(f"Organism folder already exists: {organism_folder}")
            print(f"Skipping {organism} as it has already been processed.")
            continue
        
        detailed_results, summary_df = retrieve_organism_data(
            organism=organism,
            deg_path=deg_path,
            output_path="genes/"
        )
        
        if detailed_results is not None:
            all_detailed_results[organism] = detailed_results
        
        if summary_df is not None:
            all_summaries.append(summary_df)
    
    print(f"\n{'='*50}")
    print("Script completed!")
    
    # Merge all summaries if we have multiple organisms
    if len(all_summaries) > 1:
        print("\nMerging all organism summaries...")
        merge_all_summaries("genes/")
    
    # Print final summary
    for organism, df in all_detailed_results.items():
        pdb_count = df[
            (df['pdb'].notna()) & 
            (df['pdb'] != 'Error') & 
            (df['pdb'] != 'None')
        ].shape[0]
        
        print(f"\n{organism}:")
        print(f"  Proteins: {len(df)}")
        print(f"  With downloaded PDB structures: {pdb_count}")
        print(f"  Success rate: {(pdb_count/len(df)*100):.1f}%")
    
    print(f"\n{'='*50}")
    print("All tasks completed!")

if __name__ == "__main__":
    main()