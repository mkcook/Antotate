import argparse
import tellurium as te
from equilibrator_api import ComponentContribution
from bioservices import *

class Annotate:
    def __init__(self):
        self.cc = ComponentContribution()
        self.rhea = Rhea()
    
    def parse_input_file(self, file_path):
        """Reads a text file and extracts enzyme-reaction mappings."""
        enzyme_reactions = {}
        with open(file_path, 'r') as file:
            content = file.readlines()
        
        for line in content:
            if "->" in line:
                parts = line.split(":")
                if len(parts) > 1:
                    reaction_part = parts[1].strip().replace(";", "")
                    left, right = reaction_part.split("->")
                    left_words = set(left.split())
                    right_words = set(right.split())
                    common_words = left_words.intersection(right_words)
                    enzymes = [word for word in common_words if not (word.isdigit() or word in ["+"])]
                    
                    if enzymes:
                        enzyme = enzymes[0]  # Assuming one enzyme per reaction
                        species = (left_words | right_words) - set([enzyme, "+", "->"])
                        enzyme_reactions[enzyme] = list(species)
        return enzyme_reactions
    
    def pull_tell_specs(self, file_path):
        """Pulls all floating species IDs from Antimony file."""
        with open(file_path, 'r') as file:
            antimony = file.read()
        r  = te.loada(antimony)
        return r.getFloatingSpeciesIds()
    
    def extract_chebi_ids(self, species_list):
        """Extracts ChEBI IDs using equilibrator API."""
        chebi_ids = []
        for spc in species_list:
            for i in self.cc.search_compound(spc).identifiers:
                if i.registry.namespace == 'chebi':
                    chebi_ids.append(f"{i.accession}")
                    break  # Take the first ChEBI ID found
        return chebi_ids
    
    def get_rhea_reaction(self, chebi_ids):
        """Queries Rhea for a reaction based on ChEBI IDs."""
        rh = Rhea(verbose=False)
        try:
            query = " ".join(chebi_ids)
            print(query)
            data = rh.search(query, limit=3, columns='rhea-id')
            if not data.empty:
                return data.iloc[0][0]  # Return the first matching Rhea ID
        except Exception as e:
            print(f"Error querying Rhea: {e}")
        return None
    
    def annotate_species(self, species_list, database):
        """Annotates species using the specified database."""
        annotations = []
        for spc in species_list:
            relid = None
            identity = None

            for i in self.cc.search_compound(spc).identifiers:
                if relid is None and i.registry.namespace == 'metacyc.compound':
                    relid = i.accession
                if identity is None and i.registry.namespace == database:
                    identity = i.accession
                if relid is not None and identity is not None:
                    break
            annotations.append((spc, relid, identity))
        return annotations
    
    def write_annotations(self, file_path, annotations, enzyme_annotations, database):
        """Appends annotations to the input file."""
        db_links = {
            'kegg': "http://identifiers.org/kegg/{}",
            'bigg.metabolite': "http://bigg.ucsd.edu/universal/metabolites/{}",
            'chebi': "https://www.ebi.ac.uk/chebi/searchId.do?chebiId={}",
            'hmdb': "https://hmdb.ca/metabolites/{}",
            'metacyc.compound': "https://metacyc.org/compound?orgid=META&id={}"
        }
        
        with open(file_path, 'a') as file:
            file.write('\n')  # Add a newline before appending
            for spc, relid, identity in annotations:
                file.write(f'{spc} is "{relid}";\n')
                if identity:
                    file.write(f'{spc} identity "{db_links.get(database, "http://identifiers.org/kegg/{}").format(identity)}";\n')
            for enzyme, rhea_id in enzyme_annotations.items():
                file.write(f'{enzyme} is "{enzyme}";\n')
                if rhea_id:
                    file.write(f'{enzyme} identity "https://www.rhea-db.org/rhea/{rhea_id}";\n')
                
    def annotate(self, file_path, database='kegg'):
        """Runs the annotation pipeline and returns the annotated model as text."""
        # Extract enzyme-reaction mappings
        enzyme_reactions = self.parse_input_file(file_path)
        
        # Process Antimony file
        all_species_list = self.pull_tell_specs(file_path)
        
        # Annotate species
        annotations = self.annotate_species(all_species_list, database)
        
        # Extract ChEBI IDs and query Rhea for enzyme reactions
        enzyme_annotations = {}
        for enzyme, species in enzyme_reactions.items():
            chebi_ids = self.extract_chebi_ids(species)
            rhea_id = self.get_rhea_reaction(chebi_ids) if chebi_ids else None
            enzyme_annotations[enzyme] = rhea_id
        
        # Write annotations
        self.write_annotations(file_path, annotations, enzyme_annotations, database)
        
        # Return annotated model as text
        with open(file_path, 'r') as file:
            return file.read()
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate species and enzymes in a text file.")
    parser.add_argument('file_path', help="Path to the input text file.")
    parser.add_argument(
        '--database', 
        choices=['kegg', 'bigg.metabolite', 'chebi', 'hmdb', 'metacyc.compound'], 
        default='kegg', 
        help="Annotation database to use (default: kegg)."
    )
    args = parser.parse_args()
    
    annotator = Annotate()
    annotated_text = annotator.annotate(args.file_path, args.database)
    print(annotated_text)