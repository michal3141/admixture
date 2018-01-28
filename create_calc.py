#!/usr/bin/env python

#./create_calc.py <CALC_NAME>
# The <CALC_NAME>.conf for configuration file by default

import sys
import os
import yaml
import random

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

class Individual(object):
    def __init__(self, id, sex, ethnicity):
        self.id = id
        self.sex = sex
        self.ethnicity = ethnicity

    def __str__(self):
        return 'Individual(id=%s, sex=%s, ethnicity=%s)' % (self.id, self.sex, self.ethnicity)

    def __repr__(self):
        return self.__str__()

class AdmixCalc(object):
    def __init__(self, calc_name, calc_config):
        self.calc_name = calc_name
        self.components = calc_config['components']
        self.number_of_components = calc_config['number_of_components']
        self.dataset = calc_config['dataset']
        self.admixture_params = calc_config['admixture_params']
        self.data_management = calc_config['data_management']
        self.snp_file = calc_config.get('snp_file', None)
        self.pfile = '%s.%d.P' % (self.calc_name, self.number_of_components)
        self.qfile = '%s.%d.Q' % (self.calc_name, self.number_of_components)

    def convert_eigenstrat_to_packedped(self):
        if self.data_management['convert_eigenstrat_to_packedped']:
            print 'Creating a par file: %s' % self.dataset
            with open('%s.par' % self.dataset, 'w') as par_file:
                par_file.write('genotypename: %s.geno\n' % self.dataset)
                par_file.write('snpname: %s.snp\n' % self.dataset)
                par_file.write('indivname: %s.ind\n' % self.dataset)
                par_file.write('outputformat: PACKEDPED\n')
                par_file.write('genotypeoutname: %s.bed\n' % self.dataset)
                par_file.write('snpoutname: %s.bim\n' % self.dataset)
                par_file.write('indivoutname: %s.fam\n' % self.dataset)
            print 'Converting .geno, .snp and .ind into .bed, .bim and .fam files.'
            os.system('convertf -p %s.par' % self.dataset)
        else:
            print 'Skipping coversion from eigenstrat to packedped for %s dataset.' % self.dataset

    def get_individuals(self):
        print 'Reading a %s.ind file.' % self.dataset
        self.individuals = {}
        individuals = []
        self.population_to_component = {}
        with open('%s.ind' % self.dataset, 'r') as ind_file:
            for line in ind_file:
                arr = line.strip().split()
                individ_id = arr[0]
                individ_sex = arr[1]
                individ_ethnicity = arr[2]
                individuals.append(Individual(individ_id, individ_sex, individ_ethnicity))

        print 'Creating a list of individuals from reference populations.'
        for component_name in self.components:
            component = self.components[component_name]
            for population in component:
                self.population_to_component[population] = component_name
                sample_size = component[population]
                selected_individuals = random.sample([i for i in individuals if i.ethnicity==population], sample_size)
                for individual in selected_individuals:
                    self.individuals[individual.id] = individual

        print self.individuals
        print self.population_to_component

    def filter_individuals(self):
        print 'Creating %s.fam %s.pop files with filtered individuals' % (self.calc_name, self.calc_name)
        with open('%s.fam' % self.dataset, 'r') as dataset_fam_file, open('%s.fam' % self.calc_name, 'w') as calc_fam_file, open('%s.pop' % self.calc_name, 'w') as calc_pop_file:
            for line in dataset_fam_file:
                arr = line.strip().split()
                id = arr[1]
                if id in self.individuals:
                    calc_fam_file.write(line)
                    calc_pop_file.write('%s\n' % self.population_to_component[self.individuals[id].ethnicity])

        print 'Choosing a subset of individuals from %s.bed and creating %s.bed' % (self.dataset, self.calc_name)
        os.system('plink --bfile %s --keep %s.fam --make-bed --out %s --noweb' % (self.dataset, self.calc_name, self.calc_name))

    def extract_snps(self):
        if self.snp_file:
            print 'Extracting only relevant SNPs from %s' % self.snp_file
            os.system('plink --bfile %s --extract %s --make-bed --out %s --noweb' % (self.calc_name, self.snp_file, self.calc_name))

        else:
            print 'Skipping extraction of SNPs because no SNP file was provided.'

    def run_admixture(self):
        if self.admixture_params['supervised']:
            print 'Running a supervised ADMIXTURE with K=%d components and %s.bed file.' % (self.number_of_components, self.calc_name)
            os.system('admixture --supervised %s.bed %d' % (self.calc_name, self.number_of_components))
        else:
            print 'Running an unsupervised ADMIXTURE with K=%d components and %s.bed file' % (self.number_of_components, self.calc_name)
            os.system('admixture %s.bed %d' % (self.calc_name, self.number_of_components))

    def prepare_calc_files(self):
        print 'Preparing %s.frq - a frequency file.' % self.calc_name
        os.system('plink --bfile %s --freq --out %s --noweb' % (self.calc_name, self.calc_name))

        print 'Preparing %s.alleles %s.F files' % (self.calc_name, self.calc_name)
        os.system("tail -n +2 %s.frq | awk '{ print $2,$3,$4 }' > unsorted_%s.alleles" % (self.calc_name, self.calc_name))
        os.system("paste -d' ' unsorted_%s.alleles %s | sort > merged_alleles_and_F_%s.txt" % (self.calc_name, self.pfile, self.calc_name))
        os.system("cut -d ' ' -f1,2,3 merged_alleles_and_F_%s.txt > %s.alleles" % (self.calc_name, self.calc_name))
        os.system("cut -d ' ' -f4- merged_alleles_and_F_%s.txt > %s.F" % (self.calc_name, self.calc_name))

        print 'Preparing %s.txt file' % self.calc_name
        os.system("awk '!seen[$0]++' %s.pop > %s.txt" % (self.calc_name, self.calc_name))

        print 'Preparing %s.par file' % self.calc_name
        with open('%s.par' % self.calc_name, 'w') as par_file:
            par_file.write('1d-7\n')
            par_file.write('%d\n' % self.number_of_components)
            par_file.write('genotype.txt\n')
            par_file.write('%d\n' % file_len('%s.alleles' % self.calc_name))
            par_file.write('%s.txt\n' % self.calc_name)
            par_file.write('%s.F\n' % self.calc_name)
            par_file.write('%s.alleles\n' % self.calc_name)
            par_file.write('verbose\n')
            par_file.write('genomewide\n')

    def prepare_calc_bundle(self):
        print 'Copying all the needed required files to a subdirectory: %s' % self.calc_name
        os.system('rm -rf %s' % self.calc_name)
        os.mkdir(self.calc_name)
        os.system('cp README.txt %s' % self.calc_name)
        os.system('cp standardize.r %s' % self.calc_name)
        os.system('cp DIYDodecadLinux64 %s' % self.calc_name)
        os.system('cp DIYDodecadWin.exe %s' % self.calc_name)
        os.system('cp genotype.txt %s' % self.calc_name)
        os.system('cp %s.alleles %s' % (self.calc_name, self.calc_name))
        os.system('cp %s.F %s' % (self.calc_name, self.calc_name))
        os.system('cp %s.par %s' % (self.calc_name, self.calc_name))
        os.system('cp %s.txt %s' % (self.calc_name, self.calc_name))


def main():
    calc_name = sys.argv[1]
    with open('%s.conf' % calc_name, 'r') as config_file:
        calc_config = yaml.load(config_file)
    calc = AdmixCalc(calc_name, calc_config)
    calc.convert_eigenstrat_to_packedped()
    calc.get_individuals()
    calc.filter_individuals()
    calc.extract_snps()
    calc.run_admixture()
    calc.prepare_calc_files()
    calc.prepare_calc_bundle()

if __name__ == '__main__':
    main()
