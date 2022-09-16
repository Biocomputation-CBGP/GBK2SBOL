import argparse
from cProfile import run
import os
from Bio import SeqIO
from rdflib import Graph, DCTERMS, Literal, URIRef
from sbol_generator import SBOLGenerator
from pysbolgraph.SBOL2Serialize import serialize_sboll2
from pysbolgraph.SBOL2Graph import SBOL2Graph

from identifiers import identifiers
from validate_sbol import validate_sbol
blacklist_features = ["source", "exon"]

dna_o = identifiers.roles.DNA
generator = SBOLGenerator()
gbk_map = {
    "DNA": identifiers.roles.DNA,
    "PROMOTER": identifiers.roles.promoter,
    "RBS": identifiers.roles.rbs,
    "CDS": identifiers.roles.cds,
    "TERMINATOR": identifiers.roles.terminator,
    "GENE": identifiers.roles.gene
}
existing_entities = []

def convert(gbk, out="",prefix="https://synbiohub.org/public/igem/", run_validator=True):
    graph = Graph()
    for record in SeqIO.parse(gbk, "genbank"):
        seq = record.seq
        name = os.path.basename(gbk.split(".")[0])
        subject = _build_uri(prefix,name)
        sources = [f for f in record.features if f.type == "source"]
        props = ""
        components = []
        sas = []

        for source in sources:
            for k, v in source.qualifiers.items():
                if isinstance(v, list):
                    props += " - " .join(v)
                else:
                    props += f' - {v}'
        if len(props) > 0:
            props = {DCTERMS.description: Literal(props)}
        else:
            props = {}

        for f in record.features:
            if f.type in blacklist_features:
                continue
            location = f.location
            complement = identifiers.roles.inline if location.strand == 1 else identifiers.roles.reverse
            quals = f.qualifiers
            if "label" in quals:
                e_name =  quals["label"][0]
            elif "note" in quals:
                e_name =  quals["note"][0]
            elif "gene" in quals:
                e_name =  quals["gene"][0]
            elif "primer" in quals:
                e_name = quals["primer"][0].split()[0]
            else:
                e_name = f.type

            cd_subject = _build_uri(prefix,e_name,default=f.type)
            c_subject = _build_uri(prefix,e_name,suffix="_c",default=f.type)
            s_subject = _build_uri(prefix,e_name,suffix="_annotation",default=f.type)

            c_subject = generator.build_children_uri(subject, c_subject)
            s_subject = generator.build_children_uri(subject, s_subject)
            role = _get_role(f)
            start = Literal(location.start)+1
            end = Literal(location.end)

            _add_cd(graph, cd_subject, dna_o, role)
            _add_component(graph, c_subject, cd_subject,subject)
            _add_sa(graph, s_subject, start, end, complement, c_subject,subject)

            components.append(c_subject)
            sas.append(s_subject)
        
        sequence_uri = _build_uri(prefix,name,suffix="_sequence")
        _add_sequence(graph, sequence_uri,subject, seq)
        _add_cd(graph, subject, dna_o,
                components=components, sas=sas, props=props)

    sbol = _serialise(graph)
    out_fn = os.path.join(out, gbk.split(os.path.sep)[-1].split(".")[0]+".xml")

    with open(out_fn, 'w') as o:
        o.write(sbol)

    if run_validator:
        res = validate_sbol(out_fn)
        if res["valid"]:
            print("No Errors")
        else:
            print(res["errors"])
    return out_fn


def _get_role(feature):
    descs = [feature.type]
    for q in feature.qualifiers.values():
        if isinstance(q, list):
            for e in q:
                descs += e.split()
        else:
            descs += q.split()
    for d in descs:
        d = d.upper()
        if d in gbk_map:
            return gbk_map[d]


def _build_uri(prefix,name,suffix=None,default=None):
    name = _get_name(name)
    if name.isdigit() and default is not None:
        name = default
    uri = f'{prefix}{name}{suffix if suffix is not None else ""}/1'
    if uri not in existing_entities:
        existing_entities.append(uri)
        return URIRef(uri)
    orig_uri = uri
    count = 0
    while uri in existing_entities:
        uri = f'{orig_uri[0:-2]}{str(count)}/1'
        count +=1
    existing_entities.append(uri)
    return URIRef(uri) 

def _add_sequence(graph, subject,cd, elements):
    elements = Literal(elements)
    for triple in generator.sequence(subject, elements, identifiers.objects.naseq):
        graph.add(triple)
    graph.add((cd, identifiers.predicates.sequence, subject))
    return graph


def _add_cd(graph, subject, c_type, c_role=None, components=[], sas=[], props={}):
    subject = URIRef(subject)
    for triple in generator.component_definition(subject, c_type, c_role, components, sas, props):
        graph.add(triple)
    return graph


def _add_component(graph, subject, definition,parent):
    subject = URIRef(subject)
    definition = URIRef(definition)
    for triple in generator.component(subject, definition,parent):
        graph.add(triple)
    return graph


def _add_sa(graph, subject, start, end, strand, c_subject,parent):
    subject = URIRef(subject)
    c_subject = URIRef(c_subject)
    for triple in generator.sequence_annotation(subject, start, end, strand, parent,c_subject):
        graph.add(triple)
    return graph


def _get_name(name):
    blacklist_chars = ["(",")",".","+"]
    for c in blacklist_chars:
        name = name.replace(c,"")
    return name.replace(" ", "_").replace("-", "_")


def _serialise(graph):
    pysbolG = SBOL2Graph()
    pysbolG += graph
    return serialize_sboll2(pysbolG).decode("utf-8")



def language_processor_args():
    parser = argparse.ArgumentParser(description="SBOL Enhancer Tool")
    parser.add_argument('filename', default=None, nargs='?',help="File to parse as Input")
    parser.add_argument('-f', '--output', help="The name of the output file",default="")
    parser.add_argument('-p', '--prefix', help="The prefix of translated objects.",default="https://synbiohub.org/public/igem/")
    parser.add_argument('-v', '--validator', help="Runs the tool without validation.",default=True)
    return  parser.parse_args()

if __name__ == "__main__":
    args = language_processor_args()
    convert(args.filename,args.output,args.prefix,args.validator)