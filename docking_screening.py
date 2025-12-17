import warnings
import os
import sys
import time
import subprocess
import secrets
import zipfile
import glob
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd

# Molecular modeling and bioinformatics libraries
from Bio import BiopythonWarning
from Bio.PDB import PDBParser, PDBIO, Select, PDBExceptions
from openbabel import pybel
import MDAnalysis as mda
from vina import Vina

# Suppress specific warnings
warnings.simplefilter('ignore', PDBExceptions.PDBConstructionWarning)
warnings.simplefilter('ignore', BiopythonWarning)

# ============================================================================
# LIGAND PREPARATION FUNCTIONS
# ============================================================================

def prepare_ligand_from_smiles(smiles, output_dir="temp_ligand"):
    """
    Convert SMILES to PDBQT format for docking
    
    Parameters:
        smiles: SMILES string of the ligand
        output_dir: Directory to save prepared ligand
        
    Returns:
        Path to the prepared PDBQT ligand file
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Create a temporary file with SMILES
    temp_smiles = os.path.join(output_dir, "ligand.smi")
    with open(temp_smiles, 'w') as f:
        f.write(smiles)
    
    # Use openbabel to convert SMILES to PDBQT
    try:
        mol = next(pybel.readfile("smi", temp_smiles))
        
        # Generate 3D coordinates
        mol.make3D()
        
        # Save as PDBQT
        ligand_pdbqt = os.path.join(output_dir, "ligand.pdbqt")
        mol.write("pdbqt", ligand_pdbqt, overwrite=True)
        
        print(f"Ligand prepared from SMILES and saved as: {ligand_pdbqt}")
        return ligand_pdbqt
        
    except Exception as e:
        print(f"Error preparing ligand from SMILES: {e}")
        return None

def calculate_ligand_center_of_mass(ligand_pdbqt):
    """
    Calculate the center of mass for a ligand using Pybel
    
    Parameters:
        ligand_pdbqt: Path to ligand PDBQT file
        
    Returns:
        List of [x, y, z] coordinates of center of mass
    """
    mol = next(pybel.readfile("pdbqt", ligand_pdbqt))
    atom_positions = np.array([atom.coords for atom in mol.atoms])
    center_of_mass = np.mean(atom_positions, axis=0)
    return center_of_mass.tolist()

# ============================================================================
# RECEPTOR PREPARATION FUNCTIONS
# ============================================================================

class ProteinSelect(Select):
    """BioPython Select class for selecting only protein residues"""
    def accept_residue(self, residue):
        # Only accept standard amino acid residues (protein residues)
        return residue.id[0] == ' '

def convert_to_pdb(input_file: str, output_pdb: str):
    """
    Convert various file formats to PDB using Open Babel.
    
    Parameters:
        input_file: Input file with extension: ENT, XYZ, PQR, MCIF, MMCIF, PDBQT
        output_pdb: Output PDB filename
    
    Returns:
        None
    """
    available_files = ('.ent', '.xyz', '.pqr', '.mcif', '.mmcif', '.pdbqt')
    
    # Determine the input format based on the file extension
    extension = os.path.splitext(input_file)[1].lower()
        
    if extension not in available_files:
        raise ValueError(f'Unsupported file format: {extension}. Please use an available format from {available_files}')
    
    # Use Open Babel to convert to PDB
    command = f'obabel -i{extension.lstrip(".")} "{input_file}" -opdb -O "{output_pdb}"'

    try:
        subprocess.run(command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print(f'Conversion successful: {input_file} -> {output_pdb}')
    except subprocess.CalledProcessError as e:
        print(f'Error during conversion: {e}')
        raise

def clean_and_convert_pdb_to_pdbqt(input_receptor: str, output_dir=None):
    """
    Clean the PDB file by removing non-protein residues and convert it to PDBQT format.
    
    Parameters:
        input_receptor: Path to input receptor file
        output_dir: Directory to save output files (if None, save in input directory)
    
    Returns:
        Tuple: (cleaned_pdb_path, output_pdbqt_path)
    """
    # Get the extension
    input_extension = os.path.splitext(input_receptor)[1].lower()
    
    # Determine output directory
    if output_dir is None:
        output_dir = os.path.dirname(input_receptor)
    else:
        os.makedirs(output_dir, exist_ok=True)

    # Verify which type of extension file
    if input_extension == '.pdb':
        # It's a PDB file, use it directly
        input_pdb = input_receptor
    elif input_extension in ('.ent', '.xyz', '.pqr', '.mcif', '.mmcif', '.mol2'):
        # Convert to PDB first
        base_name = os.path.splitext(os.path.basename(input_receptor))[0]
        input_pdb = os.path.join(output_dir, f"{base_name}.pdb")
        try:
            convert_to_pdb(input_receptor, input_pdb)
        except Exception as e:
            print(f"Error converting {input_receptor} to PDB: {e}")
            raise
    else:
        raise ValueError(f"Unsupported file format: {input_extension}")
    
    # Define filenames for cleaned PDB and output PDBQT
    base_name = os.path.splitext(os.path.basename(input_pdb))[0]
    cleaned_pdb = os.path.join(output_dir, f"{base_name}_cleaned.pdb")
    output_pdbqt = os.path.join(output_dir, f"{base_name}.pdbqt")
    
    # Parse the PDB file and remove non-protein residues
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', input_pdb)

    # Save the cleaned protein-only PDB
    io = PDBIO()
    io.set_structure(structure)
    io.save(cleaned_pdb, select=ProteinSelect())
    
    # Command to convert the cleaned PDB to PDBQT using Open Babel
    command = f'obabel -ipdb "{cleaned_pdb}" -opdbqt -O "{output_pdbqt}" -xr -p 7.4 --partialcharge eem'

    # Execute the command using subprocess
    try:
        print(f"Running command: {command}")
        subprocess.run(command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print(f'Receptor prepared: {output_pdbqt}')
    except subprocess.CalledProcessError as e:
        raise ValueError(f"Error during conversion: {e}")
    
    return cleaned_pdb, output_pdbqt

def calculate_center_of_mass(pdb_file):
    """
    Calculate center of mass of a protein from PDB file
    
    Parameters:
        pdb_file: Path to PDB file
        
    Returns:
        List of [x, y, z] coordinates of center of mass
    """
    parser = PDBParser(QUIET=True)
    receptor_structure = parser.get_structure('receptor', pdb_file)
    atoms = [atom for atom in receptor_structure.get_atoms()]
    center = [sum(coord) / len(atoms) for coord in zip(*[atom.coord for atom in atoms])]
    return [float(coord) for coord in center]

# ============================================================================
# GRID BOX CALCULATION FUNCTIONS
# ============================================================================

def compute_bounding_box(pdb_file):
    """
    Calculate the bounding box of a protein
    
    Parameters:
        pdb_file: Path to PDB file
        
    Returns:
        numpy array of [size_x, size_y, size_z]
    """
    u = mda.Universe(pdb_file)
    positions = u.atoms.positions
    min_coords = np.min(positions, axis=0)
    max_coords = np.max(positions, axis=0)
    size = max_coords - min_coords
    return size

def adjust_box_size(size, padding=5.0):
    """
    Increase the size of the bounding box by adding padding
    
    Parameters:
        size: Original box size [x, y, z]
        padding: Padding to add to each dimension
        
    Returns:
        Adjusted box size [x, y, z]
    """
    adjusted_size = size + padding
    return [float(dim) for dim in adjusted_size]

# ============================================================================
# LIGAND ALIGNMENT FUNCTION
# ============================================================================

def align_ligand_receptor(ligand_pdbqt, receptor_center):
    """
    Align ligand to receptor center of mass
    
    Parameters:
        ligand_pdbqt: Path to ligand PDBQT file
        receptor_center: [x, y, z] coordinates of receptor center
        
    Returns:
        None (ligand file is modified in place)
    """
    ligand_center = calculate_ligand_center_of_mass(ligand_pdbqt)
    translation_vector = np.array(receptor_center) - np.array(ligand_center)
    
    pybel_mol = next(pybel.readfile("pdbqt", ligand_pdbqt))
    for atom in pybel_mol.atoms:
        new_coords = np.add(atom.coords, translation_vector)
        atom.OBAtom.SetVector(*new_coords)
    
    pybel_mol.write("pdbqt", ligand_pdbqt, overwrite=True)
    print(f'Ligand centered to receptor: {ligand_pdbqt}')

# ============================================================================
# DOCKING FUNCTIONS - MODIFIED TO STORE ONLY FIRST POSE
# ============================================================================

def run_vina_blind_docking(receptor_pdbqt, 
                           ligand_pdbqt,
                           output_folder, 
                           center, 
                           size,
                           ligand_minimized=False,
                           exhaustiveness=32,
                           num_modes=8,
                           output_filename="output/"):
    """
    Runs blind docking using the vina Python library.
    Modified to only store the first binding pose.
    
    Parameters:
        receptor_pdbqt: Path to receptor PDBQT file
        ligand_pdbqt: Path to ligand PDBQT file
        output_folder: Directory to save results
        center: Grid center [x, y, z]
        size: Grid size [x, y, z]
        ligand_minimized: Whether ligand is already minimized
        exhaustiveness: Exhaustiveness parameter for Vina
        num_modes: Number of poses to generate (only 1 will be saved)
        output_filename: Base name for output files
        
    Returns:
        Tuple: (results_dict, error_message)
    """
    vina_results = {}
    v = Vina(sf_name='vina', cpu=-1)
    
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)

    successful_grid = False
    max_attempts = 10
    attempt = 0

    print("\nSetting up docking grid...")

    while not successful_grid and attempt < max_attempts:
        try:
            print(f"  Attempt {attempt+1}: Box size {[round(s, 1) for s in size]}")
            v.compute_vina_maps(center=center, box_size=size)
            energy = v.score()  
            successful_grid = True
        except RuntimeError as e:
            print(f"  Error: {str(e)[:100]}...")
            size = [dim + 15.0 for dim in size]
            attempt += 1
    
    if not successful_grid:
        return None, f"Failed to set the grid box after {max_attempts} attempts."

    vina_results['grid_box'] = size
    vina_results['grid_center'] = center
    vina_results['minimization'] = {}

    if not ligand_minimized:
        energy = v.score()
        print(f'  Score before minimization: {energy[0]:.3f} kcal/mol')
        vina_results['minimization']['before_min'] = energy[0]
        
        energy_minimized = v.optimize()
        print(f'  Score after minimization: {energy_minimized[0]:.3f} kcal/mol')
        vina_results['minimization']['after_min'] = energy_minimized[0]
        v.write_pose(ligand_pdbqt, overwrite=True)
    else:
        vina_results['minimization']['before_min'] = None
        vina_results['minimization']['after_min'] = None

    # Perform docking
    start_time = time.time()
    try:
        v.dock(exhaustiveness=exhaustiveness, n_poses=num_modes)
    except RuntimeError as e:
        print(f"Docking failed: {e}")
        return None, str(e)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"  Docking time: {elapsed_time:.2f} seconds")
    
    vina_results["elapsed_time"] = elapsed_time
    vina_results['results'] = []

    # Get only the first (best) pose
    try:
        energies = v.energies(n_poses=num_modes)
        print("\n  Docking completed. Best pose:")
        
        # FIX: Check if we got any results
        if energies is not None and len(energies) > 0:
            energy = energies[0]
            print(f"    Pose 1: Affinity = {energy[0]:.2f} kcal/mol")
            vina_results['results'].append({
                "pose": 1, 
                "affinity": float(energy[0]), 
                "rmsd_lb": float(energy[1]) if len(energy) > 1 else None, 
                "rmsd_ub": float(energy[2]) if len(energy) > 2 else None
            })
            
            # Store best affinity
            vina_results['best_affinity'] = float(energy[0])
        else:
            print("    No poses generated")
            vina_results['best_affinity'] = None
    except Exception as e:
        print(f"    Error getting energies: {e}")
        vina_results['best_affinity'] = None

    # Save ONLY the first pose (if available)
    os.makedirs(output_folder, exist_ok=True)
    
    # FIX: Check if best_affinity exists and is not None
    if vina_results.get('best_affinity') is not None:
        try:
            # Save the best pose
            best_pose_file = os.path.join(output_folder, "best_pose.pdbqt")
            v.write_poses(best_pose_file, n_poses=1, overwrite=True)
            vina_results['best_pose_file'] = best_pose_file
            
            # Save complex (receptor + best pose)
            complex_file = os.path.join(output_folder, "best_complex.pdbqt")
            with open(complex_file, 'w') as out_f:
                # Write receptor
                with open(receptor_pdbqt, 'r') as rec_f:
                    out_f.write(rec_f.read())
                
                # Write best pose (skip REMARK lines)
                with open(best_pose_file, 'r') as pose_f:
                    for line in pose_f:
                        if not line.startswith('REMARK'):
                            out_f.write(line)
            
            vina_results['complex_file'] = complex_file
            
            print(f"  Best pose saved: {best_pose_file}")
            print(f"  Complex saved: {complex_file}")
        except Exception as e:
            print(f"  Error saving pose files: {e}")
            vina_results['best_pose_file'] = None
            vina_results['complex_file'] = None
    else:
        vina_results['best_pose_file'] = None
        vina_results['complex_file'] = None

    print("  Results saved (no zip file created to save space)")

    return vina_results, None

# ============================================================================
# PROTEIN PROCESSING FUNCTIONS
# ============================================================================

def process_protein(pdb_file, ligand_pdbqt, organism_folder, output_base_dir, 
                    exhaustiveness=128, num_modes=1):  # Only 1 mode
    """
    Process a single protein for docking
    
    Parameters:
        pdb_file: Path to protein PDB file
        ligand_pdbqt: Path to ligand PDBQT file
        organism_folder: Organism name for folder organization
        output_base_dir: Base output directory
        exhaustiveness: Vina exhaustiveness parameter
        num_modes: Number of docking poses (only 1 will be saved)
        
    Returns:
        Dictionary with docking results or None if failed
    """
    protein_name = os.path.splitext(os.path.basename(pdb_file))[0]
    output_folder = os.path.join(output_base_dir, organism_folder, protein_name)
    
    print(f"\n{'='*60}")
    print(f"Processing protein: {protein_name}")
    print(f"{'='*60}")
    
    try:
        # Prepare receptor
        print("  Preparing receptor...")
        cleaned_pdb, receptor_pdbqt = clean_and_convert_pdb_to_pdbqt(
            pdb_file, 
            output_folder
        )
        
        if not receptor_pdbqt or not os.path.exists(receptor_pdbqt):
            print(f"  Failed to prepare receptor for {protein_name}")
            return None
        
        # Calculate box parameters
        print("  Calculating grid box...")
        center = calculate_center_of_mass(cleaned_pdb)
        size = adjust_box_size(compute_bounding_box(cleaned_pdb))
        
        # Align ligand to receptor center
        print("  Aligning ligand to receptor...")
        align_ligand_receptor(ligand_pdbqt, center)
        
        # Run docking
        print("  Running docking...")
        results, error = run_vina_blind_docking(
            receptor_pdbqt=receptor_pdbqt,
            ligand_pdbqt=ligand_pdbqt,
            output_folder=output_folder,
            center=center,
            size=size,
            ligand_minimized=False,
            exhaustiveness=exhaustiveness,
            num_modes=num_modes,
            output_filename=f"{protein_name}_results"
        )
        
        if error:
            print(f"  Docking failed: {error}")
            return None
        
        # Add protein information to results
        results['protein_name'] = protein_name
        results['pdb_file'] = os.path.basename(pdb_file)
        results['output_folder'] = output_folder
        results['receptor_pdbqt'] = receptor_pdbqt
        
        print(f"  - {protein_name} completed successfully")
        return results
        
    except Exception as e:
        print(f"  - ERROR processing {protein_name}: {str(e)[:200]}")
        return None

def create_summary_csv(results_list, output_csv):
    """
    Create a summary CSV file from docking results
    
    Parameters:
        results_list: List of docking result dictionaries
        output_csv: Path to output CSV file
        
    Returns:
        pandas DataFrame with summary or None if no results
    """
    summary_data = []
    
    for result in results_list:
        if result is None:
            continue
            
        summary_entry = {
            'Protein': result.get('protein_name', 'Unknown'),
            'PDB_File': result.get('pdb_file', 'Unknown'),
            'Best_Affinity_(kcal/mol)': result.get('best_affinity', None),
            'Grid_Center_X': result.get('grid_center', [None, None, None])[0],
            'Grid_Center_Y': result.get('grid_center', [None, None, None])[1],
            'Grid_Center_Z': result.get('grid_center', [None, None, None])[2],
            'Grid_Size_X': result.get('grid_box', [None, None, None])[0],
            'Grid_Size_Y': result.get('grid_box', [None, None, None])[1],
            'Grid_Size_Z': result.get('grid_box', [None, None, None])[2],
            'Docking_Time_(s)': result.get('elapsed_time', None),
            'Output_Folder': result.get('output_folder', ''),
            'Best_Pose_File': result.get('best_pose_file', ''),
            'Complex_File': result.get('complex_file', '')
        }
        
        # Add pose information (only pose 1)
        pose_results = result.get('results', [])
        if pose_results and len(pose_results) > 0:
            pose_result = pose_results[0]
            summary_entry['Pose_1_Affinity'] = pose_result.get('affinity', None)
            summary_entry['Pose_1_RMSD_LB'] = pose_result.get('rmsd_lb', None)
            summary_entry['Pose_1_RMSD_UB'] = pose_result.get('rmsd_ub', None)
        
        summary_data.append(summary_entry)
    
    if summary_data:
        df = pd.DataFrame(summary_data)
        df = df.sort_values('Best_Affinity_(kcal/mol)', ascending=True)
        df.to_csv(output_csv, index=False)
        print(f"\nSummary CSV saved to: {output_csv}")
        return df
    else:
        print("\nNo valid results to save to CSV")
        return None

# ============================================================================
# MAIN FUNCTION - FIXED PARAMETERS
# ============================================================================

def main():
    """Main function for high-throughput docking screening"""
    # Configuration
    DIPICOLINIC_ACID_SMILES = "C1=CC(=NC(=C1)C(=O)O)C(=O)O"
    PROTEINS_FOLDER = "genes/bacillus_subtilis_168/"
    OUTPUT_BASE_DIR = "results"
    ORGANISM_FOLDER = "bacillus_subtilis_168"
    
    # Testing parameters
    TEST_MODE = True
    MAX_TEST_PROTEINS = 20
    
    # Docking parameters
    EXHAUSTIVENESS = 32
    NUM_MODES = 8
    
    # Create output directory structure
    print("\n" + "="*60)
    print("AUTOMATED DOCKING SCREENING PIPELINE")
    print("="*60)
    print("Mode: Only storing best pose")
    print(f"Exhaustiveness: {EXHAUSTIVENESS}")
    print(f"Poses per protein: {NUM_MODES}")
    print("="*60)
    
    if not os.path.exists(OUTPUT_BASE_DIR):
        os.makedirs(OUTPUT_BASE_DIR)
        print(f"Created results directory: {OUTPUT_BASE_DIR}")
    
    organism_output_dir = os.path.join(OUTPUT_BASE_DIR, ORGANISM_FOLDER)
    os.makedirs(organism_output_dir, exist_ok=True)
    print(f"Created organism directory: {organism_output_dir}")
    
    # Prepare ligand from SMILES
    print(f"\nPreparing ligand from SMILES: {DIPICOLINIC_ACID_SMILES}")
    ligand_pdbqt = prepare_ligand_from_smiles(
        DIPICOLINIC_ACID_SMILES, 
        output_dir=os.path.join(organism_output_dir, "ligand_prep")
    )
    
    if not ligand_pdbqt or not os.path.exists(ligand_pdbqt):
        print("Failed to prepare ligand. Exiting.")
        return
    
    print(f"Ligand prepared successfully: {ligand_pdbqt}")
    
    # Get list of PDB files
    pdb_files = sorted(glob.glob(os.path.join(PROTEINS_FOLDER, "*.pdb")))
    
    if not pdb_files:
        print(f"\nNo PDB files found in {PROTEINS_FOLDER}")
        return
    
    print(f"\nFound {len(pdb_files)} PDB files in {PROTEINS_FOLDER}")
    
    # Limit to first N files for testing if in test mode
    if TEST_MODE:
        pdb_files = pdb_files[:MAX_TEST_PROTEINS]
        print(f"Running in TEST MODE - processing first {len(pdb_files)} proteins")
    
    # Process each protein
    all_results = []
    
    for i, pdb_file in enumerate(pdb_files, 1):
        print(f"\n\nProcessing protein {i}/{len(pdb_files)}")
        
        results = process_protein(
            pdb_file=pdb_file,
            ligand_pdbqt=ligand_pdbqt,
            organism_folder=ORGANISM_FOLDER,
            output_base_dir=OUTPUT_BASE_DIR,
            exhaustiveness=EXHAUSTIVENESS,
            num_modes=NUM_MODES
        )
        
        all_results.append(results)
    
    # Create summary CSV
    summary_csv = os.path.join(organism_output_dir, "docking_summary.csv")
    print(f"\n{'='*60}")
    print("CREATING SUMMARY")
    print(f"{'='*60}")
    
    df_summary = create_summary_csv(all_results, summary_csv)
    
    if df_summary is not None:
        print("\nTop 5 best docking results:")
        top_5 = df_summary[['Protein', 'Best_Affinity_(kcal/mol)']].head()
        for _, row in top_5.iterrows():
            print(f"  {row['Protein']}: {row['Best_Affinity_(kcal/mol)']:.2f} kcal/mol")
    
    print(f"\n{'='*60}")
    print("SCREENING COMPLETED!")
    print(f"{'='*60}")
    print(f"Results saved in: {organism_output_dir}")
    print(f"Summary CSV: {summary_csv}")
    
    successful = len([r for r in all_results if r is not None])
    failed = len([r for r in all_results if r is None])
    
    print(f"Total proteins processed: {len(pdb_files)}")
    print(f"Successful dockings: {successful}")
    print(f"Failed dockings: {failed}")
    
    # Calculate average affinity of successful dockings
    if successful > 0:
        affinities = [r['best_affinity'] for r in all_results if r is not None and r['best_affinity'] is not None]
        if affinities:
            avg_affinity = sum(affinities) / len(affinities)
            print(f"Average affinity: {avg_affinity:.2f} kcal/mol")
            print(f"Strongest binding: {min(affinities):.2f} kcal/mol")
            print(f"Weakest binding: {max(affinities):.2f} kcal/mol")
    
    print(f"{'='*60}")

# ============================================================================
# ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    main()