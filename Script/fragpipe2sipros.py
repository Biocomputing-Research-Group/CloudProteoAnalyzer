import getopt, sys, os
import re
import sys
# import sipros_post_module
import glob
# divide = sipros_post_module.divide
FDR_parameter = 1.0

## Glboal variables
pep_file_ext = '.pep.txt'
psm_file_ext = '.psm.txt'

pep_iden_str = '[Peptide_Identification]'
fasta_database_str = 'FASTA_Database'
pro_iden_str = '[Protein_Identification]'
testing_decoy_prefix_str = 'Testing_Decoy_Prefix'
training_decoy_prefix_str = 'Training_Decoy_Prefix'
min_peptide_per_protein_str = 'Min_Peptide_Per_Protein'
min_unique_peptide_per_protein_str = 'Min_Unique_Peptide_Per_Protein'
remove_decoy_identification_str = 'Remove_Decoy_Identification'
search_type_str = 'Search_Type'

sipros4_psmout_column_length = 20
sipros4_input = None

# defaul value
decoy_prefix = 'DECOY_'
min_peptide_per_protein = 1
min_unique_peptide_per_protein = 1
remove_decoy_identification = 'No'

class PSM:
    def __init__(self, filename, file, scan, ParentCharge, MeasuredParentMass,CalculatedParentMass,Massdiff, score, PTM_score, IdentifiedPeptide,PSM_Label,
                 Proteins,Proteinname,ProteinCount):
        self.filename=filename
        self.file = file
        self.scan = scan
        self.ParentCharge = ParentCharge
        self.MeasuredParentMass=MeasuredParentMass
        self.CalculatedParentMass=CalculatedParentMass
        self.Massdiff = Massdiff
        self.MassErrorPPM='NA'
        self.ScanType='HCD'
        self.SearchName='pFind'
        self.ScoringFunction='softmax'
        self.score = score
        self.DeltaZ='NA'
        self.DeltaP='NA'
        self.PTM_score = PTM_score
        self.IdentifiedPeptide = IdentifiedPeptide
        self.OriginalPeptide='NA'
        self.PSM_Label = PSM_Label
        self.Proteins = Proteins
        self.Proteinname=Proteinname
        self.ProteinCount=ProteinCount


class Peptide:
    def __init__(self):
        self.IdentifiedPeptide = ''
        self.ParentCharge = ''
        self.OriginalPeptide = ''
        self.ProteinNames = []
        self.ProteinCount = 0
        self.SpectralCount = 0
        self.BestScore = 0.0
        self.PSMs = []

    def add(self, oPsm):
        self.SpectralCount += 1
        if self.BestScore < oPsm.score:
            self.BestScore = oPsm.score
        self.PSMs.append('{0}_{1}'.format(oPsm.file, oPsm.scan))
        self.ScanType = 'HCD'
        self.SearchName = 'pFind'
        if oPsm.PSM_Label==True:
            self.TargetMatch='T'
        else:
            self.TargetMatch='F'

    def set(self, oPsm):
        self.IdentifiedPeptide = oPsm.IdentifiedPeptide
        self.ParentCharge = oPsm.ParentCharge
        self.OriginalPeptide = oPsm.OriginalPeptide
        self.ProteinNames = oPsm.Proteinname
        self.ProteinCount = oPsm.ProteinCount
        self.SpectralCount = 1
        self.BestScore = oPsm.score
        self.PSMs.append('{0}_{1}'.format(oPsm.file, oPsm.scan))
        self.ScanType = 'HCD'
        self.SearchName = 'pFind'
        if oPsm.PSM_Label==True:
            self.TargetMatch = 'T'
        else:
            self.TargetMatch = 'F'


def show_Fdr(psm_list, fdr_float, charge_left_given=-1, charge_right_given=-1):
    # list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
    list_sorted = sorted(psm_list, key=lambda psm: (psm.score, psm.Massdiff, psm.PTM_score), reverse=True)
    decoy = 0
    target = 0
    best_nums = [0, 0]

    psm_filtered_list = []
    cutoff_probability = 1000.0
    # without considering training label
    for oPsm in list_sorted:
        '''
        if charge_left_given != -1 and (
                oPsm.ParentCharge < charge_left_given or oPsm.ParentCharge > charge_right_given):
            continue
        '''
        if oPsm.PSM_Label:
            target += 1
        elif not oPsm.PSM_Label:
            decoy += 1
        else:
            sys.stderr.write('error 768.\n')
        (FDR_accept, FDR_value) = FDR_calculator(decoy, target)
        if (FDR_accept is True) and (FDR_value <= fdr_float):
            if (best_nums[0] + best_nums[1]) < (decoy + target):
                best_nums = [decoy, target]
                cutoff_probability = oPsm.score

    for oPsm in list_sorted:
        if charge_left_given != -1 and (
                oPsm.ParentCharge < charge_left_given or oPsm.ParentCharge > charge_right_given):
            continue
        if oPsm.score >= cutoff_probability:
            psm_filtered_list.append(oPsm)

    print("PSM level:"+str(best_nums[0])+'\t'+str(best_nums[1]))
    return psm_filtered_list


## peptide level filtering
def show_Fdr_Pep(psm_list, fdr_float, charge_left_given=-1, charge_right_given=-1):
    list_sorted = sorted(psm_list, key=lambda psm: (psm.score, psm.Massdiff, psm.PTM_score), reverse=True)

    peptide_set = set()
    decoy = 0
    target = 0
    best_nums = [0, 0]

    psm_filtered_list = []
    cutoff_probability = 1000.0
    # without considering training label
    for oPsm in list_sorted:
        '''
        if charge_left_given != -1 and (
                oPsm.ParentCharge < charge_left_given or oPsm.ParentCharge > charge_right_given):
            continue
        '''
        pep_str = oPsm.IdentifiedPeptide + '_' + str(oPsm.ParentCharge)
        if pep_str not in peptide_set:
            if oPsm.PSM_Label:
                target += 1
                peptide_set.add(pep_str)
            elif not oPsm.PSM_Label:
                decoy += 1
                peptide_set.add(pep_str)
            else:
                sys.stderr.write('error 768.\n')

        (FDR_accept, FDR_value) = FDR_calculator(decoy, target)
        if (FDR_accept is True) and (FDR_value <= fdr_float):
            if (best_nums[0] + best_nums[1]) < (decoy + target):
                best_nums = [decoy, target]
                cutoff_probability = oPsm.score

    peptide = dict()
    for oPsm in list_sorted:
        '''
        if charge_left_given != -1 and (
                oPsm.ParentCharge < charge_left_given or oPsm.ParentCharge > charge_right_given):
            continue
        '''
        pep_str=oPsm.IdentifiedPeptide+'_'+str(oPsm.ParentCharge)
        #pep_str = oPsm.IdentifiedPeptide
        if oPsm.score >= cutoff_probability:
            if pep_str in peptide:
                peptide[pep_str].add(oPsm)
            else:
                oPeptide=Peptide()
                oPeptide.set(oPsm)
                peptide[pep_str]=oPeptide

    # return set(psm_filtered_list)
    print("Peptide level:"+str(best_nums[0])+'\t'+str(best_nums[1]))
    return peptide


def read_prophet(input_file, mix_version=False):
    PSMs=[]
    C_pattern = re.compile('C\[160\]')
    M_pattern = re.compile('M\[147\]')
    clean_pattern = re.compile('[">/]')
    scan_id = 0
    charge_id = ''
    original_pep = ''
    identified_pep = ''
    protein_l = []
    iProbability = 0.0
    ntt = 0
    nmc = 0
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("<spectrum_query "):
                count=0
                split_l = line.split(' ')
                for one in split_l:
                    if one.startswith('spectrum='):
                        split_l_2 = one.split('=')
                        filename=split_l_2[-1].split('.')[0].replace('"','')+'.ms2'
                        file_id=split_l_2[-1].split('.')[0].replace('"','')
                    if one.startswith('start_scan='):
                        split_l_2 = one.split('=')
                        scan_id = int(clean_pattern.sub('', split_l_2[-1]))
                    if one.startswith('precursor_neutral_mass='):
                        split_l_2 = one.split('=')
                        MeasuredParentMass= float(clean_pattern.sub('', split_l_2[-1]))
                    if one.startswith('assumed_charge='):
                        split_l_2 = one.split('=')
                        charge_id = clean_pattern.sub('', split_l_2[-1])

                protein_l = []
                ntt = 2
                nmc = 0
            if line.startswith("<parameter name=\"ntt\""):
                split_l = line.split(' ')
                for one in split_l:
                    if one.startswith('value='):
                        split_l_2 = one.split('=')
                        ntt = int(clean_pattern.sub('', split_l_2[-1]))
            if line.startswith("<parameter name=\"nmc\""):
                split_l = line.split(' ')
                for one in split_l:
                    if one.startswith('value='):
                        split_l_2 = one.split('=')
                        nmc = int(clean_pattern.sub('', split_l_2[-1]))

            if line.startswith("<search_hit"):
                count+=1
                if count>1:
                    continue
                split_l = line.split(' ')
                for one in split_l:
                    if one.startswith('peptide='):
                        split_l_2 = one.split('=')
                        original_pep = clean_pattern.sub('', split_l_2[-1])
                        identified_pep = original_pep
                        PTM_score=identified_pep.count('~')
                    if one.startswith("protein="):
                        split_l_2 = one.split('=')
                        protein_l.append(clean_pattern.sub('', split_l_2[-1]))
                    if one.startswith("calc_neutral_pep_mass="):
                        split_l_2 = one.split('=')
                        CalculatedParentMass=float(clean_pattern.sub('', split_l_2[-1]))
                    if one.startswith("massdiff="):
                        split_l_2 = one.split('=')
                        MassDiff=float(clean_pattern.sub('', split_l_2[-1]))

            if line.startswith("<modification_info modified_peptide"):
                split_l = line.split(' ')
                for one in split_l:
                    if one.startswith('modified_peptide='):
                        split_l_2 = one.split('=')
                        identified_pep = C_pattern.sub('C', (clean_pattern.sub('', split_l_2[-1])))
                        identified_pep = M_pattern.sub('M~', (clean_pattern.sub('', split_l_2[-1])))
            if line.startswith("<alternative_protein"):
                split_l = line.split(' ')
                for one in split_l:
                    if one.startswith('protein='):
                        split_l_2 = one.split('=')
                        tmp_str = clean_pattern.sub('', split_l_2[-1])
                        if tmp_str not in protein_l:
                            protein_l.append(tmp_str)
            if line.startswith("<peptideprophet_result "):
                split_l = line.split(' ')
                for one in split_l:
                    if one.startswith('probability='):
                        split_l_2 = one.split('=')
                        iProbability = float(clean_pattern.sub('', split_l_2[-1]))
            if line.startswith("</spectrum_query>"):
                PSM_Label = False
                for p in protein_l:
                    if 'DECOY_' not in p:
                        PSM_Label = True
                        break
                PSMs.append(PSM(filename,file_id,scan_id,charge_id,MeasuredParentMass, CalculatedParentMass, MassDiff,iProbability,PTM_score,identified_pep,PSM_Label,protein_l,','.join(protein_l),len(protein_l)))
                protein_l = []


    return PSMs

# # Exit system with error message
def die(msg=None):
    if msg is not None:
        print >> sys.stderr, msg
        sys.exit(1)
# # Division error handling
def divide(x, y):
    try:
        result = x / y
    except ZeroDivisionError as detail:
        print >> sys.stderr, 'Handling run-time error:', detail
        die('Program exit!')
    else:
        return result


# # FDR calculator
def FDR_calculator(FP, TP):
    FDR_numerator = float(FP) * float(FDR_parameter)
    FDR_denominator = float(TP)
    FDR_accept = True

    if FDR_denominator == 0:
        FDR_value = 1.0
        FDR_accept = False
    else:
        FDR_value = divide(FDR_numerator, FDR_denominator)
        FDR_accept = True

    return (FDR_accept, float(FDR_value))


## remove redundant psm, only one unique spectrum kept
def re_rank(psm_list, consider_charge_bool=False):
    psm_new_list = []
    psm_dict = {}
    if consider_charge_bool:
        for oPsm in psm_list:
            sId = '{0}_{1}_{2}'.format(str(oPsm.file), str(oPsm.scan), str(oPsm.ParentCharge))
            if sId in psm_dict:
                if oPsm.score > psm_dict[sId].score:
                    psm_dict[sId] = oPsm
                elif oPsm.score == psm_dict[sId].score:
                    if abs(oPsm.Massdiff) < abs(psm_dict[sId].Massdiff):
                        psm_dict[sId] = oPsm
                    elif abs(oPsm.Massdiff) == abs(psm_dict[sId].Massdiff):
                        # calculate PTM scores
                        if oPsm.PTM_score < psm_dict[sId].PTM_score:
                            psm_dict[sId] = oPsm
                        elif oPsm.PTM_score == psm_dict[sId].PTM_score:
                            if oPsm.IdentifiedPeptide.upper() < psm_dict[sId].IdentifiedPeptide.upper():
                                psm_dict[sId] = oPsm
                            elif oPsm.IdentifiedPeptide.upper() == psm_dict[sId].IdentifiedPeptide.upper():
                                psm_dict[sId].add_protein(oPsm.protein_list)

            else:
                psm_dict[sId] = oPsm
    else:
        for oPsm in psm_list:
            sId = '{0}_{1}'.format(str(oPsm.file), str(oPsm.scan))
            if sId in psm_dict:
                if oPsm.score > psm_dict[sId].score:
                    psm_dict[sId] = oPsm
                elif oPsm.score == psm_dict[sId].score:
                    if abs(oPsm.Massdiff) < abs(psm_dict[sId].Massdiff):
                        psm_dict[sId] = oPsm
                    elif abs(oPsm.Massdiff) == abs(psm_dict[sId].Massdiff):
                        # calculate PTM scores
                        if oPsm.PTM_score < psm_dict[sId].PTM_score:
                            psm_dict[sId] = oPsm
                        elif oPsm.PTM_score == psm_dict[sId].PTM_score:
                            if oPsm.IdentifiedPeptide.upper() < psm_dict[sId].IdentifiedPeptide.upper():
                                psm_dict[sId] = oPsm


            else:
                psm_dict[sId] = oPsm

    for _key, value in psm_dict.items():
        psm_new_list.append(value)

    return psm_new_list



if __name__ == "__main__":
    '''
    file_str=sys.argv[1]
    file_str_prefix=file_str.split('.')[0]
    PSMs=readpFind(file_str)
    '''
    file_folder=sys.argv[1]
    files=glob.glob(file_folder+'*.pep.xml')
    file_str_prefix=file_folder.split('/')[-2]
    fdr=float(sys.argv[2])
    PSMs=[]
    for file in files:
        PSMs.extend(read_prophet(file))
    rank_list = re_rank(PSMs)
    filter_list = show_Fdr(rank_list, fdr)
    filter_pep_list = show_Fdr_Pep(rank_list, fdr)
    

    with open(file_str_prefix+'.psm.txt', 'w') as f:
        psm_out_list = ['Filename',  # 0
                    'ScanNumber',  # 1
                    'ParentCharge',  # 2
                    'MeasuredParentMass',  # 3
                    'CalculatedParentMass',  # 4
                    'MassErrorDa',  # 5 CalculatedParentMass - MeasuredParentMass
                    'MassErrorPPM',  # 6 MassErrorDa / CalculatedParentMass
                    'ScanType',  # 7
                    'SearchName',  # 8
                    'ScoringFunction',  # 9
                    'Score',  # 10
                    'DeltaZ',  # 11 the difference score between the rank 1 and 2
                    'DeltaP',  # 12
                    'IdentifiedPeptide',  # 13
                    'OriginalPeptide',  # 14
                    'ProteinNames',  # 15
                    'ProteinCount',  # 16
                    'TargetMatch']  # 17
        f.write('\t'.join(psm_out_list) + '\n')

            
        for psm in filter_list:
            TargetMatch = 'F'
            if psm.PSM_Label == True:
                TargetMatch = 'T'
            f.write(str(psm.filename)+ '\t' + str(psm.scan) + '\t' + str(psm.ParentCharge)+'\t'+str(psm.MeasuredParentMass)+'\t'+str(psm.CalculatedParentMass)+'\t'+str(psm.Massdiff)+'\t'+str(psm.MassErrorPPM)+'\t'+str(psm.ScanType)+'\t'+str(psm.SearchName)+'\t'+str(psm.ScoringFunction)+'\t'+str(-psm.score)+'\t'+str(psm.DeltaZ)+'\t'+str(psm.DeltaP)+'\t'+ str(psm.IdentifiedPeptide) + '\t' +str(psm.OriginalPeptide)+'\t'+str(psm.Proteinname)+'\t'+ str(psm.ProteinCount)+'\t'+TargetMatch + '\n')
    

    with open(file_str_prefix+'.pep.txt', 'w') as f:
        pep_out_list = ['IdentifiedPeptide',    #0
                    'ParentCharge',         #1
                    'OriginalPeptide',      #2
                    'ProteinNames',         #3
                    'ProteinCount',         #4
                    'TargetMatch',          #5
                    'SpectralCount',        #6 number of PSMs matched to this peptide
                    'BestScore',            #7 the highest score of those PSMs
                    'PSMs',                 #8 a list of PSMs matched to this peptide. Use{Filename[ScanNumber],Filename[ScanNumber]} format
                    'ScanType',             #9
                    'SearchName']           #10
        f.write('\t'.join(pep_out_list) + '\n')
        for key, pep in filter_pep_list.items():
            f.write(pep.IdentifiedPeptide+'\t'+str(pep.ParentCharge)+'\t'+pep.OriginalPeptide+'\t'+'{'+pep.ProteinNames+'}'+'\t'+str(pep.ProteinCount)+'\t'+pep.TargetMatch+'\t'+str(pep.SpectralCount)+'\t'+str(-pep.BestScore)+'\t'+','.join(pep.PSMs)+'\t'+pep.ScanType+'\t'+pep.SearchName+ '\n')



    print('Done')