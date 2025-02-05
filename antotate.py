import argparse
import tellurium as te
import csv
from equilibrator_api import ComponentContribution
from bioservices import Rhea

class Annotate:
    def __init__(self):
        self.cc = ComponentContribution()
        self.rhea = Rhea()
    
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
    
    def extract_chebi_ids(self, species_list):
        """Extracts ChEBI IDs using equilibrator API."""
        chebi_ids = []
        for spc in species_list:
            for i in self.cc.search_compound(spc).identifiers:
                if i.registry.namespace == 'chebi':
                    chebi_ids.append(i.accession)
                    break  # Take the first ChEBI ID found
        return chebi_ids
    
    def get_rhea_reaction(self, chebi_ids):
        """Queries Rhea for a reaction based on ChEBI IDs."""
        try:
            reaction = self.rhea.get_reaction_by_chebi(chebi_ids)
            if reaction:
                return reaction[0]  # Return the first matching Rhea ID
        except Exception as e:
            print(f"Error querying Rhea: {e}")
        return None
    
    def calculate_confidence_score(self, species):
        """Calculates a confidence score for species annotation."""
        search_results = self.cc.ccache.search(species, page=1, page_size=25)
        result = self.cc.search_compound(species)
        return next((score for compound, score in search_results if compound.id == result.id), None)
    
    def annotate_species(self, species_list, databases):
        """Annotates species using the specified databases."""
        annotations = []
        confidence_scores = []
        
        for spc in species_list:
            relid = None
            identities = {}
            score = self.calculate_confidence_score(spc)

            for i in self.cc.search_compound(spc).identifiers:
                if relid is None and i.registry.namespace == 'metacyc.compound':
                    relid = i.accession
                if i.registry.namespace in databases:
                    identities[i.registry.namespace] = i.accession
            
            annotations.append((spc, relid, identities))
            confidence_scores.append((spc, relid, identities, score))
        
        return annotations, confidence_scores
    
    def write_annotations(self, file_path, annotations, enzyme_annotations, databases):
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
            for spc, relid, identities in annotations:
                file.write(f'{spc} is "{relid}";\n')
                for db, identity in identities.items():
                    file.write(f'{spc} identity "{db_links.get(db, "http://identifiers.org/kegg/{}").format(identity)}";\n')
            for enzyme, rhea_id in enzyme_annotations.items():
                file.write(f'{enzyme} is "{enzyme}";\n')
                if rhea_id:
                    file.write(f'{enzyme} identity "https://www.rhea-db.org/rhea/{rhea_id}";\n')
    
    def write_confidence_metrics(self, confidence_scores):
        """Writes confidence scores to a CSV file."""
        with open("confidence_metrics.csv", 'w', newline='') as csvfile:
            fieldnames = ['Species', 'relid', 'Identities', 'Confidence Score']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for species, relid, identities, score in confidence_scores:
                writer.writerow({'Given ID': species, 'Display Name': relid, 'Annotated Identities': ', '.join(f'{db}:{id}' for db, id in identities.items()), 'Confidence Score': score})
    
    def annotate(self, file_path, databases):
        """Runs the annotation pipeline and returns the annotated model as text."""
        # Process input file
        enzyme_list = self.parse_input_file(file_path)
        
        # Process Antimony file
        all_species_list = self.pull_tell_specs(file_path)
        
        # Remove enzymes from species list
        species_list = list(set(all_species_list) - set(enzyme_list))
        
        # Annotate species
        annotations, confidence_scores = self.annotate_species(species_list, databases)
        
        # Extract ChEBI IDs and query Rhea for enzyme reactions
        chebi_ids = self.extract_chebi_ids(species_list)
        rhea_id = self.get_rhea_reaction(chebi_ids) if chebi_ids else None
        enzyme_annotations = {enzyme: rhea_id for enzyme in enzyme_list if rhea_id}
        
        # Write annotations
        self.write_annotations(file_path, annotations, enzyme_annotations, databases)
        
        # Write confidence metrics
        self.write_confidence_metrics(confidence_scores)
        
        # Return annotated model as text
        with open(file_path, 'r') as file:
            return file.read()
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate species and enzymes in a text file.")
    parser.add_argument('file_path', help="Path to the input text file.")
    parser.add_argument('--databases', nargs='+', choices=['kegg', 'bigg.metabolite', 'chebi', 'hmdb', 'metacyc.compound'], default=['kegg'], help="Annotation databases to use (default: kegg).")
    args = parser.parse_args()
    
    annotator = Annotate()
    annotated_text = annotator.annotate(args.file_path, args.databases)
    print(annotated_text)