#!/usr/bin/env python

import sys
import argparse
import random


_SEX = {
    '0': 'UNKNOWN',
    '1': 'MALE',
    '2': 'FEMALE'
}

_KAR_SEX = {
    '0': None,
    '1': 'XY',
    '2': 'XX'
}
_VARIABLE_FIELDS = ['Relationship', 'Siblings', 'Second Order', 'Third Order', 'Other Comments']

_FNAME_TEMPLATE = 'ALL.chr{}_GRCh38_sites.20170504.vcf.gz'

_DISORDERS = [
    {'id': 'OMIM:126200', 'name': 'MULTIPLE SCLEROSIS, SUSCEPTIBILITY TO; MS', 'source': 'OMIM'},
    {'id': 'OMIM:266600', 'name': 'INFLAMMATORY BOWEL DISEASE (CROHN DISEASE) 1; IBD1', 'source': 'OMIM'},
    {'id': 'OMIM:605803', 'name': 'DERMATITIS, ATOPIC, 2; ATOD2 ', 'source': 'OMIM'},
    {'id': 'OMIM:177900', 'name': 'PSORIASIS 1, SUSCEPTIBILITY TO; PSORS1', 'source': 'OMIM'},
    {'id': 'OMIM:125853', 'name': 'DIABETES MELLITUS, NONINSULIN-DEPENDENT; NIDDM ', 'source': 'OMIM'},
    {'id': 'OMIM:209850', 'name': 'AUTISM', 'source': 'OMIM'},
    {'id': 'OMIM:106300', 'name': 'SPONDYLOARTHROPATHY, SUSCEPTIBILITY TO, 1; SPDA1', 'source': 'OMIM'},
    {'id': 'OMIM:310200', 'name': 'MUSCULAR DYSTROPHY, DUCHENNE TYPE; DMD', 'source': 'OMIM'},
    {'id': 'OMIM:606232', 'name': 'PHELAN-MCDERMID SYNDROME; PHMDS', 'source': 'OMIM'}
]

_PHENOTYPES = [
    {'id': 'HP:0001324', 'name': 'Muscle weakness', 'source': 'HPO', 'status': 'OBSERVED'},
    {'id': 'HP:0100280', 'name': 'Crohn\'s disease', 'source': 'HPO', 'status': 'OBSERVED'},
    {'id': 'HP:0001047', 'name': 'Atopic dermatitis', 'source': 'HPO', 'status': 'OBSERVED'},
    {'id': 'HP:0003765', 'name': 'Psoriasiform dermatitis', 'source': 'HPO', 'status': 'OBSERVED'},
    {'id': 'HP:0005978', 'name': 'Type II diabetes mellitus', 'source': 'HPO', 'status': 'OBSERVED'},
    {'id': 'HP:0000717', 'name': 'Autism', 'source': 'HPO', 'status': 'OBSERVED'},
    {'id': 'HP:0001369', 'name': 'Arthritis', 'source': 'HPO', 'status': 'OBSERVED'},
    {'id': 'HP:0003560', 'name': 'Muscular dystrophy', 'source': 'HPO', 'status': 'OBSERVED'},
    {'id': 'HP:0000708', 'name': 'Behavioral abnormality', 'source': 'HPO', 'status': 'OBSERVED'}
]


    # Ankylosing Spondylitis:
#     id: DOID:7147
#     name: ankylosing spondylitis
#     source: Human Disease Ontology



def to_camel_case(text):
    components = text.lower().replace('_', ' ').split(' ')
    return components[0].lower() + ''.join(x.title() for x in components[1:])


def create_variable_sets(header):
    text = []
    text.append('variableSets:')
    text.append('{}- id: relation'.format(' '*2))
    text.append('{}name: relation'.format(' '*4))
    text.append('{}entities:'.format(' '*4))
    text.append('{}- INDIVIDUAL'.format(' '*6))
    text.append('{}variables:'.format(' '*4))
    for field in header:
        if field not in _VARIABLE_FIELDS:
            continue
        text.append('{}- id: {}'.format(' '*6, to_camel_case(field)))
        text.append('{}name: {}'.format(' '*8, to_camel_case(field)))
        text.append('{}type: STRING'.format(' '*8))
    return '\n'.join(text)


def create_individuals(ind_info):
    text = []
    text.append('individuals:')
    for i, ind in enumerate(ind_info):
        text.append('{}- id: {}'.format(' '*2, ind['Individual ID']))
        text.append('{}name: {}'.format(' '*4, ind['Individual ID']))
        if ind['Paternal ID'] != '0':
            text.append('{}father:'.format(' '*4))
            text.append('{}id: {}'.format(' '*6, ind['Paternal ID']))
        if ind['Maternal ID'] != '0':
            text.append('{}mother:'.format(' '*4))
            text.append('{}id: {}'.format(' '*6, ind['Maternal ID']))
        text.append('{}sex: {}'.format(' '*4, _SEX[ind['Gender']]))
        text.append('{}karyotypicSex: {}'.format(' '*4, _KAR_SEX[ind['Gender']]))
        text.append('{}population:'.format(' '*4))
        text.append('{}name: {}'.format(' '*6, ind['Population']))

        if i % 5 == 0:
            index = random.randrange(len(_DISORDERS))

            phenotype = _PHENOTYPES[index]
            text.append('{}phenotypes:'.format(' ' * 4))
            text.append('{}- id: {}'.format(' ' * 6, phenotype['id']))
            text.append('{}name: {}'.format(' ' * 8, phenotype['name']))
            text.append('{}source: {}'.format(' ' * 8, phenotype['source']))
            text.append('{}STATUS: {}'.format(' ' * 8, phenotype['status']))

            disorder = _DISORDERS[index]
            text.append('{}disorders:'.format(' ' * 4))
            text.append('{}- id: {}'.format(' ' * 6, disorder['id']))
            text.append('{}name: {}'.format(' ' * 8, disorder['name']))
            text.append('{}source: {}'.format(' ' * 8, disorder['source']))

        text.append('{}annotationSets:'.format(' '*4))
        text.append('{}- id: relation'.format(' '*6))
        text.append('{}name: relation'.format(' '*8))
        text.append('{}variableSetId: relation'.format(' '*8, ind['Population']))
        text.append('{}annotations:'.format(' '*8))
        for field in ind.keys():
            if field not in _VARIABLE_FIELDS:
                continue
            text.append('{}{}: {}'.format(' '*10, to_camel_case(field), ind[field]))
        text.append('{}samples:'.format(' '*4))
        text.append('{}- id: {}'.format(' '*6, ind['Individual ID']))
    return '\n'.join(text)


def create_samples(ind_info):
    text = []
    text.append('samples:')
    for ind in ind_info:
        text.append('{}- id: {}'.format(' '*2, ind['Individual ID']))
        text.append('{}individualId: {}'.format(' '*4, ind['Individual ID']))
    return '\n'.join(text)


def create_families(ind_info):
    families = {}
    for ind in ind_info:
        if ind['Family ID'] != ind['Individual ID']:
            families.setdefault(ind['Family ID'], []).append(ind['Individual ID'])
            for member in [ind['Individual ID'], ind['Paternal ID'], ind['Maternal ID']]:
                if member != '0':
                    families[ind['Family ID']].append(member)

    for family in families:
        families[family] = list(set(families[family]))

    text = []
    text.append('families:')
    for family in families:
        text.append('{}- id: {}'.format(' '*2, family))
        text.append('{}name: {}'.format(' '*4, family))
        text.append('{}members:'.format(' '*4))
        for member in families[family]:
            text.append('{}- id: {}'.format(' '*6, member))
    return '\n'.join(text)


def create_files():
    text = []
    text.append('files:')
    for chrom in list(range(1, 23)) + ['X', 'Y']:
        text.append('{}- name: {}'.format(' ' * 2, _FNAME_TEMPLATE.format(chrom)))
        text.append('{}path: {}'.format(' '*4, 'data'))
    return '\n'.join(text)


def _setup_argparse():
    desc = 'This script creates automatically all Python RestClients files'
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('ped_file', help='Pedigree file path')
    parser.add_argument('outfile', help='Output file path')
    args = parser.parse_args()
    return args


def main():
    random.seed(42)

    args = _setup_argparse()

    ped_fhand = open(args.ped_file, 'r')
    yml_fhand = open(args.outfile, 'w')

    header = ped_fhand.readline().strip().split('\t')
    ind_info = [{k: v for k, v in zip(header, line.strip().split('\t'))} for line in ped_fhand]

    yml_fhand.write('---\n')
    yml_fhand.write('# WARNING: AUTOGENERATED CONTENT\n')
    yml_fhand.write('id: 1000g\n')
    yml_fhand.write('name: 1000 Genomes phase 3\n')
    yml_fhand.write('description: The 1000 Genomes Project\n')
    yml_fhand.write(create_variable_sets(header) + '\n')
    yml_fhand.write(create_individuals(ind_info) + '\n')
    yml_fhand.write(create_families(ind_info) + '\n')
    yml_fhand.write(create_files() + '\n')

    ped_fhand.close()
    yml_fhand.close()


if __name__ == '__main__':
    sys.exit(main())