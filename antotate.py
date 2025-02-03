import argparse
import tellurium as te
from equilibrator_api import ComponentContribution

class Antotate:
    def __init__(self):
        self.cc = ComponentContribution()
    
    def parse_input_file(self, file_path):
        """Reads a text file and extracts unique species."""
        with open(file_path, 'r') as file:
            content = file.readlines()
        
        species_set = set()
        for line in content:
            if "->" in line:
                parts = line.split(":")
                if len(parts) > 1:
                    reaction_part = parts[1].strip().replace(";", "")
                    left, right = reaction_part.split("->")
                    left_words = set(left.split())
                    right_words = set(right.split())
                    common_words = left_words.intersection(right_words)
                    species_set.update(word for word in common_words if not (word.isdigit() or word in ["+"]))
        return list(species_set)
    
    def pull_tell_specs(self, file_path):
        """Pulls all floating species IDs from Antimony file."""
        with open(file_path, 'r') as file:
            antimony = file.read()
        r  = te.loada(antimony)
        return r.getFloatingSpeciesIds()
    
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
    
    def write_annotations(self, file_path, annotations, database):
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
                    file.write(f'{spc} identity "{db_links.get(database, "http://identifiers.org/kegg/{}" ).format(identity)}";\n')
                
    
    def annotate(self, file_path, database='kegg'):
        """Runs the annotation pipeline and returns the annotated model as text."""
        # Process input file
        enzyme_list = self.parse_input_file(file_path)
        
        # Process Antimony file
        all_species_list = self.pull_tell_specs(file_path)
        
        # Remove enzymes from species list
        species_list = list(set(all_species_list) - set(enzyme_list))
        
        # Annotate species
        annotations = self.annotate_species(species_list, database)
        
        # Write annotations
        self.write_annotations(file_path, annotations, database)
        
        # Return annotated model as text
        with open(file_path, 'r') as file:
            return file.read()
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate species in a text file.")
    parser.add_argument('file_path', help="Path to the input text file.")
    parser.add_argument(
        '--database', 
        choices=['kegg', 'bigg.metabolite', 'chebi', 'hmdb', 'metacyc.compound'], 
        default='kegg', 
        help="Annotation database to use (default: kegg)."
    )
    args = parser.parse_args()
    
    annotator = Antotate()
    annotated_text = annotator.annotate(args.file_path, args.database)
    print(annotated_text)