import os
from resources import *
import functions
debug_label=False

def generateSlurmOpenmp(path:str, params:dict, session_dict:dict)->None:
    with open(path,"w") as f:
        f.write("#!/bin/bash\n")
        if debug_label:
            f.write("#SBATCH --partition=debug\n") ################################################
        else:
            f.write("#SBATCH --partition=normal\n")
        f.write("#SBATCH --exclusive\n")
        f.write("#SBATCH --output="+session_dict['user_email']+"_%J_stdout.txt\n")
        f.write("#SBATCH --error="+session_dict['user_email']+"_%J_stderr.txt\n")
        if debug_label:
            f.write("#SBATCH --time=00:30:00\n")
        else:
            f.write("#SBATCH --time="+str(params["hour"])+":"+str(params["mins"])+":00\n")
        f.write("#SBATCH --job-name="+session_dict['user_email']+"\n")
        f.write("module load python/anaconda2-4.2.0\n")
        if functions.readUploadFile(os.path.join(user_data_path,session_dict['user_email']),"multiple")=="1":
            f.write("SiprosEnsembleOMP -w " +"raw_data" + " -c generated_sipros_config.cfg -o result\n")           
        else:
            f.write("SiprosEnsembleOMP -f " +params["ms file"] + " -c generated_sipros_config.cfg -o result\n")
        f.write("TabFile=$(python2 ~/Scripts/sipros_psm_tabulating.py -i result/ -c generated_sipros_config.cfg -o result/)\n")
        f.write("python2 ~/Scripts/sipros_ensemble_filtering.py -i ${TabFile} -c generated_sipros_config.cfg -o result/\n")
        f.write("python2 ~/Scripts/sipros_peptides_assembling.py -w result/ -c generated_sipros_config.cfg\n")
        files=functions.readUploadFile(os.path.join(user_data_path,session_dict['user_email']),"files").strip().split("~")
        for tmp in files:
            f.write("ProRata.linux -w "+tmp+"/ -c ProRataConfig.xml  -i result/*.pro2psm.txt\n")
    f.close()

def generateSlrumMPI(path:str, path2:str, params:dict, session_dict:dict)->None:
    with open(path,"w") as f:
        f.write("#!/bin/bash\n")
        if debug_label:
            f.write("#SBATCH --partition=debug\n") ################################################
        else:
            f.write("#SBATCH --partition=normal\n")
        f.write("#SBATCH --exclusive\n")
        if debug_label:
            f.write("#SBATCH --nodes=2\n")
            f.write("#SBATCH --ntasks=20\n")
            f.write("#SBATCH --ntasks-per-node=10\n")
        else:
            f.write("#SBATCH --nodes=2\n")
            f.write("#SBATCH --ntasks=20\n")
            f.write("#SBATCH --ntasks-per-node=10\n")
        f.write("#SBATCH --output="+session_dict['user_email']+"_%J_stdout.txt\n")
        f.write("#SBATCH --error="+session_dict['user_email']+"_%J_stderr.txt\n")
        if debug_label:
            f.write("#SBATCH --time=00:30:00\n")
        else:
            f.write("#SBATCH --time=47:59:00\n")
        f.write("#SBATCH --job-name="+session_dict['user_email']+"\n")
        f.write("#SBATCH --chdir=/scratch/jiancheng/"+session_dict['whole_path'].split("/")[-2]+"+"+session_dict['whole_path'].split("/")[-1]+"\n")
        f.write("module load python/anaconda2-4.2.0\n")
        f.write("source ~/miniconda3/etc/profile.d/conda.sh\n")
        f.write("conda activate mpi\n")
        f.write("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CONDA_PREFIX}/lib\n")
        f.write("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/jiancheng/condaEnsembleBIN\n")
        f.write("#export OMP_NUM_THREADS=10\n")
        if functions.readUploadFile(os.path.join(user_data_path,session_dict['user_email']),"multiple")=="1":
            f.write("mpirun /home/jiancheng/condaEnsembleBIN/SiprosEnsembleMPI -w " +"raw_data" + " -c generated_sipros_config.cfg -o result\n")           
        else:
            f.write("mpirun /home/jiancheng/condaEnsembleBIN/SiprosEnsembleMPI -f " +params["ms file"] + " -c generated_sipros_config.cfg -o result\n")
        f.write("TabFile=$(python2 ~/Scripts/sipros_psm_tabulating.py -i result/ -c generated_sipros_config.cfg -o result/)\n")
        f.write("python2 ~/Scripts/sipros_ensemble_filtering.py -i ${TabFile} -c generated_sipros_config.cfg -o result/\n")
        f.write("python2 ~/Scripts/sipros_peptides_assembling.py -w result/ -c generated_sipros_config.cfg\n")
        if params["search type"]=="SIP":
            f.write("python2 ~/Scripts/ClusterSip.py -w result/ -c generated_sipros_config.cfg\n")
    f.close()
    with open(path2,"w") as f:
        f.write("#!/bin/bash\n")
        if debug_label:
            f.write("#SBATCH --partition=debug\n") ################################################
        else:
            f.write("#SBATCH --partition=normal\n")
        f.write("#SBATCH --exclusive\n")
        f.write("#SBATCH --output="+session_dict['user_email']+"_%J_stdout.txt\n")
        f.write("#SBATCH --error="+session_dict['user_email']+"_%J_stderr.txt\n")
        if debug_label:
            f.write("#SBATCH --time=00:30:00\n")
        else:
            f.write("#SBATCH --time="+str(params["hour"])+":"+str(params["mins"])+":00\n")
        f.write("#SBATCH --job-name="+session_dict['user_email']+"\n")
        # f.write("module load python/anaconda2-4.2.0\n")
        # f.write("TabFile=$(python2 ~/Scripts/sipros_psm_tabulating.py -i result/ -c generated_sipros_config.cfg -o result/)\n")
        # f.write("python2 ~/Scripts/sipros_ensemble_filtering.py -i ${TabFile} -c generated_sipros_config.cfg -o result/\n")
        # f.write("python2 ~/Scripts/sipros_peptides_assembling.py -w result/ -c generated_sipros_config.cfg\n")
        # if params["search type"]=="SIP":
        #     f.write("python2 ~/Scripts/ClusterSip.py -w result/ -c generated_sipros_config.cfg\n")
        #     f.write("ProRata.linux -w result/ -c ProRataConfig.xml -i result/*.pro2psm.cluster.txt")
        # else:
        #     f.write("ProRata.linux -w result/ -c ProRataConfig.xml -i result/*.pro2psm.txt")
        files=functions.readUploadFile(os.path.join(user_data_path,session_dict['user_email']),"files").strip().split("~")
        for tmp in files:
            f.write("ProRata.linux -w "+tmp+"/ -c ProRataConfig.xml  -i result/*.pro2psm.txt\n")
        
#########################################
    f.close()


def generateSiprosConfigFile(config_dict:dict, whole_path:str):
    FASTA_Database = "newdatabase.fasta"
    with open(os.path.join(whole_path,"generated_sipros_config.cfg"),"w") as f:
        f.write("###########################################################\n##### Parameters for peptide identification by Sipros #####\n###########################################################\n[Peptide_Identification]\n\n")
        f.write("Search_Name = "+config_dict["Search_Name"]+"\n")
        f.write("Search_Type = "+config_dict["Search_Type"]+'\n')
        f.write("FASTA_Database = "+FASTA_Database+"\n")
        f.write("Fragmentation_Method = "+config_dict["Fragmentation_Method"]+"\n")
        f.write("Parent_Mass_Windows = "+config_dict["Parent_Mass_Windows"]+"\n")
        f.write("Search_Mass_Tolerance_Parent_Ion = "+config_dict["Search_Mass_Tolerance_Parent_Ion"]+"\n")
        f.write("Mass_Tolerance_Fragment_Ions = "+config_dict["Mass_Tolerance_Fragment_Ions"]+"\n")
        f.write("Minimum_Peptide_Length = "+config_dict["Minimum_Peptide_Length"]+"\n")
        f.write("Cleave_After_Residues = "+ config_dict["Cleave_After_Residues"]+"\n")
        f.write("Maximum_Peptide_Length = "+config_dict["Maximum_Peptide_Length"]+"\n")
        f.write("Cleave_Before_Residues = "+config_dict["Cleave_Before_Residues"]+"\n")
        f.write("Maximum_Missed_Cleavages = "+config_dict["Maximum_Missed_Cleavages"]+"\n")
        f.write("Try_First_Methionine = "+config_dict["Try_First_Methionine"]+"\n")
        f.write("Max_PTM_Count = "+config_dict["Max_PTM_Count"]+"\n")
        for i in range(len(config_dict["PTM_Values"])):
            f.write(config_dict["PTM_Values"][i]+"\n")
        f.write("##### Elemental composition of amino acid residues #####\n")
        f.write("Element_List	=	C,	H,	O,	N,	P,	S\n")
        f.write("Residue{Nterm}	=	0,	1,	0,	0,	0,	0	# N terminus\n")
        f.write("Residue{Cterm}	=	0,	1,	1,	0,	0,	0	# C terminus\n")
        f.write("Residue{J}	=	6,	11,	1,	1,	0,	0	# J is I or L\n")
        f.write("Residue{I}	=	6,	11,	1,	1,	0,	0   \n")
        f.write("Residue{L}	=	6,	11,	1,	1,	0,	0	\n")
        f.write("Residue{A}	=	3,	5,	1,	1,	0,	0	\n")
        f.write("Residue{S}	=	3,	5,	2,	1,	0,	0	\n")
        f.write("Residue{G}	=	2,	3,	1,	1,	0,	0	\n")
        f.write("Residue{V}	=	5,	9,	1,	1,	0,	0	\n")
        f.write("Residue{E}	=	5,	7,	3,	1,	0,	0	\n")
        f.write("Residue{K}	=	6,	12,	1,	2,	0,	0	\n")
        f.write("Residue{T}	=	4,	7,	2,	1,	0,	0	\n")
        f.write("Residue{D}	=	4,	5,	3,	1,	0,	0	\n")
        f.write("Residue{R}	=	6,	12,	1,	4,	0,	0	\n")
        f.write("Residue{P}	=	5,	7,	1,	1,	0,	0	\n")
        f.write("Residue{N}	=	4,	6,	2,	2,	0,	0	\n")
        f.write("Residue{F}	=	9,	9,	1,	1,	0,	0	\n")
        f.write("Residue{Q}	=	5,	8,	2,	2,	0,	0	\n")
        f.write("Residue{Y}	=	9,	9,	2,	1,	0,	0	\n")
        f.write("Residue{M}	=	5,	9,	1,	1,	0,	1	\n")
        f.write("Residue{H}	=	6,	7,	1,	3,	0,	0	\n")
        f.write("Residue{C}	=	5,	8,	2,	2,	0,	1	# Blocked Cys by IAA\n")
        f.write("#Residue{C}	=	3,	5,	1,	1,	0,	1	# Natural Cys \n")
        f.write("Residue{W}	=	11,	10,	1,	2,	0,	0	\n")
        for one_ptm in config_dict["PTM_Values"]:
            if one_ptm.find("~")!=-1:
                f.write("Residue{~}	=	0,	0,	1,	0,	0,	0,	# Oxidation or Hydroxylation\n")
            elif one_ptm.find("!")!=-1:
                f.write("Residue{!}	=	0,	-1,	1,	-1,	0,	0,	# Deamidation or Citrullination if happens at Arg\n")
            elif one_ptm.find("@")!=-1:
                f.write("Residue{@}	=	0,	1,	3,	0,	1,	0,	# Phosphorylation\n")
                f.write("Residue{>}	=	0,	1,	3,	0,	1,	0,	# Phosphorylation\n")
                f.write("Residue{<}	=	0,	1,	3,	0,	1,	0,	# Phosphorylation\n")
            elif one_ptm.find(">to1")!=-1:
                # f.write("Residue{>}	=	0,	1,	3,	0,	1,	0,	# Phosphorylation\n")
                f.write("Residue{1}	=	0,	0,	0,	0,	0,	0,	# Phosphorylation with losing HPO3\n")
            elif one_ptm.find("<to2")!=-1:
                # f.write("Residue{<}	=	0,	1,	3,	0,	1,	0,	# Phosphorylation\n")
                f.write("Residue{2}	=	0,	-2,	-1,	0,	0,	0,	# Phosphorylation with losing HPO3 and H2O\n")
            elif one_ptm.find("%")!=-1:
                f.write("Residue{%}	=	2,	2,	1,	0,	0,	0,	# Acetylation\n")
            elif one_ptm.find("^")!=-1:
                f.write("Residue{^}	=	1,	2,	0,	0,	0,	0,	# Mono-methylation\n")
            elif one_ptm.find("&")!=-1:
                f.write("Residue{&}	=	2,	4,	0,	0,	0,	0,	# Di-methylation\n")
            elif one_ptm.find("*")!=-1:
                f.write("Residue{*}	=	3,	6,	0,	0,	0,	0,	# Tri-methylation\n")
            elif one_ptm.find("(")!=-1:
                f.write("Residue{(}	=	0,	-1,	1,	1,	0,	0,	# S-Nitrosylation\n")
            elif one_ptm.find(")")!=-1:
                f.write("Residue{)}	=	0,	-1,	2,	1,	0,	0,	# Nitration\n")
            elif one_ptm.find("/")!=-1:
                f.write("Residue{/}	=	2,	3,	1,	1,	0,	0,	# IAA blocking\n")
            elif one_ptm("$")!=-1:
                f.write("Residue{$}	=	1,	2,	0,	0,	0,	1,	# beta-methylthiolation\n")
        f.write("##### Isotopic distribution of elements #####\n")
        f.write("# Carbon\n")
        f.write("Element_Masses{C} 	=	12.000000,	13.003355\n")
        f.write("Element_Percent{C} 	=	0.9893,		0.0107\n")
        f.write("# Hydrogen\n")
        f.write("Element_Masses{H} 	=	1.007825,	2.014102\n")
        f.write("Element_Percent{H} 	=	0.999885,	0.000115\n")
        f.write("# Oxygen\n")
        f.write("Element_Masses{O} 	=	15.994915,	16.999132,	17.999160\n")
        f.write("Element_Percent{O} 	=	0.99757,	0.00038,	0.00205\n")
        f.write("# Nitrogen\n")
        f.write("Element_Masses{N} 	=	14.003074,	15.000109\n")
        f.write("Element_Percent{N} 	=	0.99632,	0.00368\n")
        f.write("# Phosphorus \n")
        f.write("Element_Masses{P} 	=	30.973762\n")
        f.write("Element_Percent{P} 	=	1.0\n")
        f.write("# Sulfur\n")
        f.write("Element_Masses{S} 	=	31.972071,	32.971459,	33.967867,	34.967867, 	35.967081\n")
        f.write("Element_Percent{S} 	=	0.9493,		0.0076,		0.0429,		0.0000, 	0.0002\n")    
        f.write("###########################################################\n##### Parameters for protein identification by Sipros #####\n###########################################################\n[Protein_Identification]\n")
        f.write("Training_Decoy_Prefix = Rev1_\n")
        f.write("Testing_Decoy_Prefix = Rev_\n")
        # f.write("Reserved_Decoy_Prefix = Rev_2_\n")
        f.write("FDR_Filtering = "+config_dict["FDR_Filtering"]+"\n")
        f.write("FDR_Threshold = "+config_dict["FDR_Threshold"]+"\n")
        f.write("Min_Peptide_Per_Protein = "+config_dict["Min_Peptide_Per_Protein"]+"\n")
        f.write("Min_Unique_Peptide_Per_Protein = "+config_dict["Min_Unique_Peptide_Per_Protein"]+"\n")
        f.write("Filter_Mass_Tolerance_Parent_Ion = "+config_dict["Filter_Mass_Tolerance_Parent_Ion"]+"\n")
        f.write("Filter_Mass_Tolerance_Parent_Ion_Unit = "+config_dict["Filter_Mass_Tolerance_Parent_Ion_Unit"]+"\n")
        if config_dict["Search_Type"].find("SIP")!=-1:
            f.write("###########################################################\n##### Parameters for stable isotope probing by Sipros #####\n###########################################################\n[Stable_Isotope_Probing]\n")
            f.write("SIP_Element = "+config_dict["SIP_Element"]+"\n")
            f.write("SIP_Element_Isotope = "+str(config_dict["SIP_Element_Isotope"])+"\n")
            f.write("Maximum_Enrichment_Level = "+str(config_dict["Maximum_Enrichment_Level"])+"%\n")
            f.write("Minimum_Enrichment_Level = "+str(config_dict["Minimum_Enrichment_Level"])+"%\n")
            f.write("Enrichment_Level_Increment = "+str(config_dict["Enrichment_Level_Increment"])+"%\n")
            f.write("Clustering_Threshold = "+str(config_dict["Clustering_Threshold"])+"\n")
            f.write("Min_PSM_Per_Isotopic_Cluster = "+str(config_dict["Min_PSM_Per_Isotopic_Cluster"])+"\n")
        f.close()

def generateProrataConfigFile(config_dict:dict, whole_path:str):
    FASTA_Database = "newdatabase.fasta"
    with open(os.path.join(whole_path,"ProRataConfig.xml"),"w") as f:
        f.write("<?xml version = \"1.0\" ?>\n<CONFIG>\n\t<!-- Parameters for SIC-Forma -->\n\t<SIC_EXTRACTION>\n\t\t<RETENTION_TIME_INTERVAL>\n")
        f.write("\t\t\t<MINUTES_BEFORE_MS2>"+config_dict["MINUTES_BEFORE_MS2"]+"</MINUTES_BEFORE_MS2>\n")
        f.write("\t\t\t<MINUTES_AFTER_MS2>"+config_dict["MINUTES_AFTER_MS2"]+"</MINUTES_AFTER_MS2>\n")
        f.write("\t\t\t<MINUTES_BETWEEN_DUPLICATE_MS2>"+config_dict["MINUTES_BETWEEN_DUPLICATE_MS2"]+"</MINUTES_BETWEEN_DUPLICATE_MS2>\n")
        f.write("\t\t</RETENTION_TIME_INTERVAL>\n\t\t<MASS_TO_CHARGE_INTERVAL>\n")
        f.write("\t\t\t<PLUS_MZ_ERROR>"+config_dict["PLUS_MZ_ERROR"]+"</PLUS_MZ_ERROR>\n")
        f.write("\t\t\t<MINUS_MZ_ERROR>"+config_dict["MINUS_MZ_ERROR"]+"</MINUS_MZ_ERROR>\n")
        f.write("\t\t\t<ISOTOPIC_ENVELOP_CUTOFF>"+config_dict["ISOTOPIC_ENVELOP_CUTOFF"]+"</ISOTOPIC_ENVELOP_CUTOFF>\n")
        f.write("\t\t</MASS_TO_CHARGE_INTERVAL>\n\t\t<!-- ATOM ISOTOPIC COMPOSITION -->\n\t\t<ATOM_ISOTOPIC_COMPOSITION>\n")
        f.write("\t\t\t<!-- Data taken from http://physics.nist.gov/PhysRefData/Compositions/index.html -->\n")
        f.write("\t\t\t<C>\n")
        f.write("\t\t\t\t<MASS_DA>	12.000000,	13.003355	</MASS_DA>\n")
        f.write("\t\t\t\t<NATURAL>	0.9893,		0.0107		</NATURAL>\n")
        f.write("\t\t\t\t<ENRICHED>	0.02,		0.98		</ENRICHED>\n")
        f.write("\t\t\t</C>\n")
        f.write("\t\t\t<H>\n")
        f.write("\t\t\t\t<MASS_DA>	1.007825,	2.014102	</MASS_DA>\n")
        f.write("\t\t\t\t<NATURAL>	0.999885,	0.000115	</NATURAL>\n")
        f.write("\t\t\t\t<ENRICHED>	0.02,		0.98		</ENRICHED>\n")
        f.write("\t\t\t</H>\n")
        f.write("\t\t\t<O>\n")
        f.write("\t\t\t\t<MASS_DA>	15.994915,	16.999132,	17.999160	</MASS_DA>\n")
        f.write("\t\t\t\t<NATURAL>	0.99757,	0.00038,	0.00205		</NATURAL>\n")
        f.write("\t\t\t\t<ENRICHED>	0.02,		0.0,		0.98		</ENRICHED>\n")
        f.write("\t\t\t</O>\n")
        f.write("\t\t\t<N>\n")
        f.write("\t\t\t\t<MASS_DA>	14.003074,	15.000109	</MASS_DA>\n")
        f.write("\t\t\t\t<NATURAL>	0.99632,	0.00368		</NATURAL>\n")
        f.write("\t\t\t\t<ENRICHED>	0.02,		0.98		</ENRICHED>\n")
        f.write("\t\t\t</N>\n")
        f.write("\t\t\t<P>\n")
        f.write("\t\t\t\t<MASS_DA>	30.973762	</MASS_DA>\n")
        f.write("\t\t\t\t<NATURAL>	1.0		</NATURAL>\n")
        f.write("\t\t\t\t<ENRICHED>	1.0		</ENRICHED>\n")
        f.write("\t\t\t</P>\n")
        f.write("\t\t\t<S>\n")
        f.write("\t\t\t\t<MASS_DA>	31.972071,	32.971459,	33.967867,	35.967081	</MASS_DA>\n")
        f.write("\t\t\t\t<NATURAL>	0.9493,		0.0076,		0.0429,		0.0002		</NATURAL>\n")
        f.write("\t\t\t\t<ENRICHED>	0.02,		0.0,		0.98,		0.0		</ENRICHED>\n")
        f.write("\t\t\t</S>\n")
        f.write("\t\t</ATOM_ISOTOPIC_COMPOSITION>\n")
        f.write("\t\t<RESIDUE_ATOMIC_COMPOSITION>\n")
        f.write("\t\t\t<ISOTOPOLOGUE name = \"label_free\" >\n")
        f.write("\t\t\t\t<!--	Name	C	H	O	N	P	S	C*	H*	O*	N*	P*	S*	-->\n")
        f.write("\t\t\t\t<R>	NTerm,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n")
        f.write("\t\t\t\t<R>	CTerm,	0,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n")
        f.write("\t\t\t\t<R>	L,	6,	11,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n")
        f.write("\t\t\t\t<R>	A,	3,	5,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n")
        f.write("\t\t\t\t<R>	S,	3,	5,	2,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n")
        f.write("\t\t\t\t<R>	G,	2,	3,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n")
        f.write("\t\t\t\t<R>	V,	5,	9,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n")
        f.write("\t\t\t\t<R>	E,	5,	7,	3,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n")
        f.write("\t\t\t\t<R>	K,	6,	12,	1,	2,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n")
        f.write("\t\t\t\t<R>	I,	6,	11,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n")
        f.write("\t\t\t\t<R>	T,	4,	7,	2,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n")
        f.write("\t\t\t\t<R>	D,	4,	5,	3,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n")
        f.write("\t\t\t\t<R>	R,	6,	12,	1,	4,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n")
        f.write("\t\t\t\t<R>	P,	5,	7,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n")
        f.write("\t\t\t\t<R>	N,	4,	6,	2,	2,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n")
        f.write("\t\t\t\t<R>	F,	9,	9,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n")
        f.write("\t\t\t\t<R>	Q,	5,	8,	2,	2,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n")
        f.write("\t\t\t\t<R>	Y,	9,	9,	2,	1,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n")
        f.write("\t\t\t\t<R>	M,	5,	9,	1,	1,	0,	1,	0,	0,	0,	0,	0,	0	</R>\n")
        f.write("\t\t\t\t<R>	H,	6,	7,	1,	3,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n")
        f.write("\t\t\t\t<R>	C,	3,	5,	1,	1,	0,	1,	0,	0,	0,	0,	0,	0	</R>\n")
        f.write("\t\t\t\t<R>	W,	11,	10,	1,	2,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n")
        for one_ptm in config_dict["PTM_Values"]:
            f.write(ptm_iso_free_prorata[one_ptm])
        f.write("\t\t\t</ISOTOPOLOGUE>\n")

        if config_dict["Search_Type"].find("SIP")!=-1:
            if config_dict["SIP_Element"]=="C":
                f.write(C_iso_prorata)
                for one_ptm in config_dict["PTM_Values"]:
                    f.write(ptm_iso_c_prorata[one_ptm])
                f.write("			</ISOTOPOLOGUE>\n")
            elif config_dict["SIP_Element"]=="N":
                f.write(N_iso_prorata)
                for one_ptm in config_dict["PTM_Values"]:
                    f.write(ptm_iso_n_prorata[one_ptm])
                f.write("			</ISOTOPOLOGUE>\n")
        
        f.write("\t\t</RESIDUE_ATOMIC_COMPOSITION>\n\t</SIC_EXTRACTION>\n")
        f.write("\t<PEPTIDE_QUANTIFICATION>\n\t\t<PEAK_DETECTION>\n")
        f.write("\t\t<CHROMATOGRAM_SMOOTHING>\n")
        f.write("\t\t\t\t<ORDER>"+config_dict["ORDER"] +"</ORDER>\n")
        f.write("\t\t\t\t<WINDOW_SIZE>"+config_dict["WINDOW_SIZE"]+"</WINDOW_SIZE>\n")
        f.write("\t\t\t</CHROMATOGRAM_SMOOTHING>\n")
        f.write("\t\t\t<PEAK_SHIFT>\n")
        f.write("\t\t\t\t<LEFT>"+config_dict["LEFT"]+"</LEFT>\n")
        f.write("\t\t\t\t<RIGHT>"+config_dict["RIGHT"]+"</RIGHT>\n")
        f.write("\t\t\t</PEAK_SHIFT>\n")
        f.write("\t\t</PEAK_DETECTION>\n")
        if config_dict["Search_Type"].find("SIP")!=-1:
            f.write("		<ABUNDANCE_RATIO>\n			<NUMERATOR_ISOTOPOLOGUE>	label_free	</NUMERATOR_ISOTOPOLOGUE>\n			<DENOMINATOR_ISOTOPOLOGUE>	reference	</DENOMINATOR_ISOTOPOLOGUE>\n		</ABUNDANCE_RATIO>\n")
            f.write("		<LOG2_RATIO>\n")
            f.write("			<MINIMUM>	"+config_dict["LOGp_MINIMUM"]+"	</MINIMUM>\n")
            f.write("			<MAXIMUM>	"+config_dict["LOGp_MAXIMUM"]+"	</MAXIMUM>\n")
            f.write("		</LOG2_RATIO>\n		<LOG2_SNR_CUTOFF>	"+config_dict["LOGp_SNR_CUTOFF"]+"	</LOG2_SNR_CUTOFF>\n")
        f.write("\t\t<REMOVE_AMBIGUOUS_PEPTIDES>"+config_dict["REMOVE_AMBIGUOUS_PEPTIDES"]+"</REMOVE_AMBIGUOUS_PEPTIDES>\n")
        f.write("\t</PEPTIDE_QUANTIFICATION>\n")
        if config_dict["Search_Type"].find("SIP")!=-1:
            f.write("        	<PROTEIN_QUANTIFICATION>\n")
            f.write("		<FASTA_FILE>		"+FASTA_Database+"	</FASTA_FILE>\n")
            f.write("		<MIN_PEPTIDE_NUMBER>	"+config_dict["MIN_PEPTIDE_NUMBER"]+"	</MIN_PEPTIDE_NUMBER>\n")
            f.write("		<MAX_CI_WIDTH>		"+config_dict["MAX_CI_WIDTH"]+"	</MAX_CI_WIDTH>\n")
            f.write("		<MAX_LOG2_SNR>		"+config_dict["MAX_LOG2_SNR"]+"	</MAX_LOG2_SNR>\n")
            f.write("		<LOG2_RATIO>\n")
            f.write("			<MINIMUM>	"+config_dict["LOGT_MINIMUM"]+"	</MINIMUM>\n")
            f.write("			<MAXIMUM>	"+config_dict["LOGT_MAXIMUM"]+"	</MAXIMUM>\n")
            f.write("		</LOG2_RATIO>\n")
            f.write("		<LOG2_RATIO_DISCRETIZATION>	"+config_dict["LOG2_RATIO_DISCRETIZATION"]+"	</LOG2_RATIO_DISCRETIZATION>\n")
            f.write("		<STANDARD_DEVIATION>\n			<SLOPE>	-0.288	</SLOPE>\n			<INTERCEPT> 1.305	</INTERCEPT>\n		</STANDARD_DEVIATION>\n		<MEAN>\n			<SLOPE>	1.2	</SLOPE>\n			<INTERCEPT>	0	</INTERCEPT>\n		</MEAN>\n		<SMOOTHING_PROBABILITY_SPACE>	0.15    </SMOOTHING_PROBABILITY_SPACE>\n	</PROTEIN_QUANTIFICATION>\n")
        f.write("</CONFIG>")
    f.close()    
