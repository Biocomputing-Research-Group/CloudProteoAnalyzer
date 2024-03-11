from distutils.command.config import config
from multiprocessing.sharedctypes import Value
from nis import match
from pickle import TRUE
import streamlit as st #pip install streamlit
import os
import time
import shutil
import time #pip install times
import io
import google.oauth2.credentials #pip install --upgrade google-auth
import google_auth_oauthlib.flow #pip install google-auth-oauthlib
import googleapiclient.discovery as gdiscovery #pip install google-api-python-client
import googleapiclient as gapi
import google.auth.transport.requests #pip install google-auth
import inspect #pip install inspect-it
import google_auth_oauthlib as gauth
import requests #pip install requests
import functions
from GetAlphaPeptMS2 import *
import logging
import genconfigfiles
import pandas as pd
from resources import *
import base64
from pathlib import Path
from PIL import Image

st.set_page_config(
    page_title="CloudProteoAnalyzer",
    page_icon="SP",
    layout="wide",
    menu_items={
        'Get Help': None,
        'Report a bug': None,
        'About': "# The CloudProteoAnalyzer"
    })
# ptm_list = ["PTM{~} = M","PTM{~} = MPK","PTM{!} = NQR","PTM{@} = STYHD","PTM{>to1} = STYHD","PTM{<to2} = ST","PTM{%} = K","PTM{^} = KRED","PTM{&} = KR","PTM{*} = K","PTM{(} = C","PTM{)} = Y","PTM{/} = C","PTM{$} = D"]


def Downloadfromdrive(session_dict:dict,sipros_path:str,config_dict:dict,page_drive:st):
    flow = st.session_state.oauth_flow
    start_index = 0
    stop_index = 0
    err_label = ""
    files=""

    for i in range(len(config_dict["folderlink"])):
        if config_dict["folderlink"][i:i+len('folders/')] == 'folders/':
            start_index = i+len('folders/')
        if config_dict["folderlink"][i:i+len('?usp')] == '?usp':
            stop_index = i
    folder_id = config_dict["folderlink"][start_index:stop_index]

    drive = gdiscovery.build('drive', 'v3', credentials=flow.credentials)
    page_token = None
    found_folder = False
    all_files=[]
    all_files_name = []
    while True:
        try:
            response = drive.files().list(q= "mimeType='application/vnd.google-apps.folder'",
                                        spaces='drive',
                                        fields='nextPageToken, files(id, name)',
                                        pageToken=page_token).execute()
        except:
            st.write("Please give the permission to download from your google drive!")
            return
        for file in response.get('files', []):
            if file.get('id') == folder_id:
                found_folder = True
                page_token = response.get('nextPageToken', None)
        if page_token is None or found_folder == True:
            break
    if not found_folder:
        # st.error("Your uploading folder is not in your google drive.\n Please upload the new folder.")
        err_label="no file"
        shutil.rmtree(os.path.join(sipros_path,session_dict['user_email'],"0"))
        return False, err_label
    while True:
        response = drive.files().list(q= "'" + folder_id + "' in parents",
                                              spaces='drive',
                                              fields='nextPageToken, files(id, name)',
                                              pageToken=page_token).execute()
        for file in response.get('files', []):
            all_files.append(file.get('id'))
            all_files_name.append(file.get('name'))
            page_token = response.get('nextPageToken', None)
        if page_token is None:
            break
            
    if not functions.Checkfiles(all_files_name):
        err_label="miss file"
        shutil.rmtree(os.path.join(sipros_path,session_dict['user_email'],"0"))
        return False, err_label
    count_raw = 0
    page_drive.write('Please wait for the file to finish uploading! Do not cleck any button.')
    for one_file in range(len(all_files)):
        request = drive.files().get_media(fileId=all_files[one_file])
        fh = io.FileIO(session_dict['whole_path']+"/"+all_files_name[one_file], 'wb')
        downloader = gapi.http.MediaIoBaseDownload(fh, request)
        done = False
        while done is False:
            status, done = downloader.next_chunk()
        if all_files_name[one_file].find(".raw") != -1:
            count_raw+=1
            if all_files_name[one_file].find(".RAW") != -1:
                os.system("mv "+session_dict['whole_path']+"/"+all_files_name[one_file]+" "+session_dict['whole_path']+"/"+all_files_name[one_file].replace(".RAW",".raw"))
            os.system("mono "+raxport_path+"/Raxport.exe -i "+session_dict['whole_path']+" -o "+session_dict['whole_path'])
            files+=all_files_name[one_file][:-4]+"~"
            os.remove(session_dict['whole_path'] +"/"+all_files_name[one_file])
        # elif all_files_name[one_file].find(".fa") != -1:
        #     functions.Reversedatabase(session_dict['whole_path'] +"/"+all_files_name[one_file],session_dict['whole_path']+"/newdatabase.fasta",config_dict)
        #     os.remove(session_dict['whole_path']+"/"+all_files_name[one_file])
    if count_raw>1:
        mul_raw_data = "1"
    else:
        mul_raw_data = "0"
    page_drive.write("Uploading Finish!")
    functions.writeUploadFile(os.path.join(sipros_path,session_dict['user_email']),("",""),files[:-1])
    # print(files)
    return mul_raw_data,err_label

def paramsTable(form:st, config_dict:dict):
    label_search = form.selectbox(label="Search type:(Regular: label free quantification.)",options=["Regular"])
    config_dict["Version_selection"]=form.selectbox(label="Version_Selection:\nThere are two basic versions of the database-searching: one for running on a single machine and another for running with MPI on a cluster.",options=["single machine", "Cluster"])
    
    pept_expander=form.expander("Parameters for peptide identification")
    config_dict["Fragmentation_Method"]  = pept_expander.selectbox(label='Fragmentation_Method:',options=["CID","HCD"])#string dropdown
    config_dict["Parent_Mass_Windows"]  = pept_expander.text_input(label='Parent_Mass_Windows:\nParent mass tolerance for database searching\nMass Windows to be open around parent ion.\nExamples: a center window: "0", a center window and an offset window: "-1,0", etc',value="-1,0,1,2,3")#int slider? Not sure
    config_dict["Search_Mass_Tolerance_Parent_Ion"]   = str(pept_expander.number_input(label='Search_Mass_Tolerance_Parent_Ion\n\
                                                                     Recommend 1 Da for both High-res MS2 and Low-res MS2 to improve identification results\n\
                                                                     Peptides with large errors for parent ions can be filtered out using the parameter Filter_Mass_Tolerance_Parent_Ion\n',min_value=0.01, max_value=1.0,step=0.01,value=0.05))#int slider 0.05-0.5
    config_dict["Mass_Tolerance_Fragment_Ions"] =str(pept_expander.number_input(label="Mass_Tolerance_Fragment_Ions:\ne.g. \"0.05\" for High-res MS2 and \"0.5\" for Low-res MS2.",min_value=0.01,max_value=1.0,step=0.01,value=0.05))
    config_dict["Minimum_Peptide_Length"]  = str(int(pept_expander.number_input(label='Minimum_Peptide_Length:(1-15)',min_value=1,max_value=15,step=1,value=7)))#int slider 1-15
    config_dict["Maximum_Peptide_Length"]  = str(int(pept_expander.number_input(label='Maximum_Peptide_Length:(20-100)',min_value=20,max_value=100,step=1,value=60)))#int slider 20-100
    enzyme_type=pept_expander.selectbox(label='Enzyme_type:',options=["Trypsin","Trypsin/P","Lys_C","Lys_N","Arg_C","Asp_N","CNBr","Glu_C","PepsinA","Chymotrypsin"])
    config_dict["Cleave_After_Residues"]   = enzyme_dict[enzyme_type][0]#string textbox                
    config_dict["Cleave_Before_Residues"]   = enzyme_dict[enzyme_type][1]#string textbox
    config_dict["Maximum_Missed_Cleavages"]   = str(int(pept_expander.number_input(label='Maximum_Missed_Cleavages:(1-5)\nMaximum number of missed cleavages in a peptide',min_value=1,max_value=5,step=1,value=3)))#int slider 1-5
    config_dict["Try_First_Methionine"]   = pept_expander.selectbox(label='Try_First_Methionine:\nTry removing the first methionine in a protein',options=["True","False"])#bool dropdown
    config_dict["Max_PTM_Count"]   = str(int(pept_expander.number_input(label='Max_PTM_Count:(1-5)\nconsidered for a peptide',min_value=1,max_value=5,step=1,value=3)))#int slider 1-5
    config_dict["PTM_Values"] = pept_expander.multiselect(label="PTMs:",options=ptm_list,default=["PTM{~} = M"]) 

    ###
    prot_expander=form.expander("Parameters for protein identification")
    config_dict["FDR_Filtering"] = prot_expander.selectbox(label='FDR_Filtering:\nLevel of FDR filtering',options=["PSM","Peptide"])#string dropdown
    config_dict["FDR_Threshold"] = str(prot_expander.number_input(label='FDR_Threshold:(0.01-0.5)\nFDR threshold for filtering peptide identifications',min_value=0.01,max_value=0.05,step=0.01,value=0.01))#int slider 
    config_dict["Min_Peptide_Per_Protein"] = str(int(prot_expander.number_input(label='Min_Peptide_Per_Protein:(1-6):\nMinimum number of peptides per protein',min_value=1,max_value=6,step=1,value=1)))#int slider
    config_dict["Min_Unique_Peptide_Per_Protein"] = str(int(prot_expander.number_input(label='Min_Unique_Peptide_Per_Protein:(1-6)\nMinimum number of unique peptides per protein',min_value=1,max_value=6,step=1,value=1)))#int slider
    config_dict["Filter_Mass_Tolerance_Parent_Ion"] = str(float(prot_expander.number_input(label='Filter_Mass_Tolerance_Parent_Ion:\ne.g. "0.05" for High-res MS1 and "3" for Low-res MS1, default unit is Da',min_value=0.001,max_value=50.0,step=0.1,value=0.05)))
    
    # prot_expander.selectbox(label='Filter_Mass_Tolerance_Parent_Ion:\ne.g. "0.05" for High-res MS1 and "3" for Low-res MS1, default unit is Da',options=["0.05","3"])#int slider
    config_dict["Filter_Mass_Tolerance_Parent_Ion_Unit"] = prot_expander.selectbox(label='Filter_Mass_Tolerance_Parent_Ion_Unit:\nAtomic mass error unit',options=["Da","PPM"])#string dropdown
             
    if label_search == "SIP":
        config_dict["Search_Type"] = "SIP"
        sip_expander=form.expander("Parameters for stable isotope probing")
        config_dict["SIP_Element"] = sip_expander.selectbox(label="SIP_Element:\nThe partially labeled element.\n",options=["C","N"])
        if config_dict["SIP_Element"]=="C":
            config_dict["SIP_Element_Isotope"] = sip_expander.selectbox(label="SIP_Element_Isotope:\nIntegar mass of the isotope of the element to be searched",options=[12,13])
        # elif config_dict["SIP_Element"]=="H":
        #     config_dict["SIP_Element_Isotope"]= sip_expander.selectbox(label="SIP_Element_Isotope:\nIntegar mass of the isotope of the element to be searched",options=[1,2])
        # elif config_dict["SIP_Element"]=="O":
        #     config_dict["SIP_Element_Isotope"]= sip_expander.selectbox(label="SIP_Element_Isotope:\nIntegar mass of the isotope of the element to be searched",options=[15,16,17])
        elif config_dict["SIP_Element"]=="N":
            config_dict["SIP_Element_Isotope"]= sip_expander.selectbox(label="SIP_Element_Isotope:\nIntegar mass of the isotope of the element to be searched",options=[14,15])
        # elif config_dict["SIP_Element"]=="P":
        #     config_dict["SIP_Element_Isotope"]= sip_expander.selectbox(label="SIP_Element_Isotope:\nIntegar mass of the isotope of the element to be searched",options=[30])
        # elif config_dict["SIP_Element"]=="S":
        #     config_dict["SIP_Element_Isotope"]= sip_expander.selectbox(label="SIP_Element_Isotope:\nIntegar mass of the isotope of the element to be searched",options=[31,32,33,34,35])
        config_dict["Maximum_Enrichment_Level"] = int(sip_expander.number_input(label="Maximum_Enrichment_Level(%):\nThe maximum enrichment levels to be searched",min_value=0,max_value=100,step=1,value=100))
        config_dict["Minimum_Enrichment_Level"] = sip_expander.number_input(label="Minimum_Enrichment_Level(%):\nThe minimum enrichment levels to be searched",min_value=0,max_value=config_dict["Maximum_Enrichment_Level"]-1,step=1,value=0)
        config_dict["Enrichment_Level_Increment"] = sip_expander.number_input(label="Enrichment_Level_Increment(%):\nThe increment of the enrichment levels to be searched",min_value=0,max_value=50,step=1,value=1)
        config_dict["Clustering_Threshold"] = sip_expander.number_input(label="Clustering_Threshold(%):\nEnrichment level threshold for merging peptide clusters of a protein",min_value=0,max_value=config_dict["Maximum_Enrichment_Level"]-1,step=1,value=20)
        config_dict["Min_PSM_Per_Isotopic_Cluster"] = sip_expander.number_input(label="Min_PSM_Per_Isotopic_Cluster:\n",min_value=1,max_value=50,step=1,value=1)
        config_dict["Search_Name"]=config_dict["SIP_Element"]+"_"+str(config_dict["SIP_Element_Isotope"])
    else:
        config_dict["Search_Type"]="Regular"
        config_dict["Search_Name"]="SE"

    ###
    quant_expander=form.expander("Parameters for protein quantification")
    rtint_container=quant_expander.container()
    rtint_container.write("Retiotion time interval for SIC extraction:")
    config_dict["MINUTES_BEFORE_MS2"] =str(rtint_container.number_input(label="MINUTES_BEFORE_MS2:",min_value=0.5,max_value=3.0,step=0.5,value=2.0))
    config_dict["MINUTES_AFTER_MS2"] =str(rtint_container.number_input(label="MINUTES_AFTER_MS2:",min_value=0.5,max_value=3.0,step=0.5,value=2.0))
    config_dict["MINUTES_BETWEEN_DUPLICATE_MS2"] =str(rtint_container.number_input(label="MINUTES_BETWEEN_DUPLICATE_MS2:",min_value=0.5,max_value=3.0,step=0.5,value=2.0))

    ####
    mz_container=quant_expander.container()
    mz_container.write("")
    mz_container.write("Mass to charge interval:")
    config_dict["PLUS_MZ_ERROR"] =str(mz_container.number_input(label="PLUS_MZ_ERROR:",min_value=0.01,max_value=1.0,step=0.01,value=0.5))
    config_dict["MINUS_MZ_ERROR"] =str(mz_container.number_input(label="MINUS_MZ_ERROR:",min_value=0.01,max_value=1.0,step=0.01,value=0.5))
    config_dict["ISOTOPIC_ENVELOP_CUTOFF"] =str(mz_container.number_input(label="ISOTOPIC_ENVELOP_CUTOFF:",min_value=0.01,max_value=1.0,step=0.01,value=0.1))

    ###
    smooth_container=quant_expander.container()
    smooth_container.write("")
    smooth_container.write("Chromatogram smoothing during peak detection:")
    config_dict["ORDER"] =str(int(smooth_container.selectbox(label="ORDER:",options=[2,3])))
    config_dict["WINDOW_SIZE"] =str(int(smooth_container.number_input(label="WINDOW_SIZE:",min_value=5,max_value=10,step=1,value=7)))

    ###
    shift_container=quant_expander.container()
    shift_container.write("")
    shift_container.write("Peak shift during peak detection:")
    config_dict["LEFT"] =str(shift_container.number_input(label="PEAK_SHIFT_LEFT:",min_value=0,max_value=3,step=1,value=0))
    config_dict["RIGHT"] =str(shift_container.number_input(label="PEAK_SHIFT_RIGHT:",min_value=0,max_value=3,step=1,value=0))
    quant_expander.write("")

    if label_search == "SIP":
        Peptquant_container=quant_expander.container()
        Peptquant_container.write("")
        Peptquant_container.write("Log2 ratio for peptide quantification:")
        config_dict["LOGp_MINIMUM"]=str(Peptquant_container.number_input(label="MINIMUM:",min_value=-10,max_value=10,step=1,value=-10))
        config_dict["LOGp_MAXIMUM"]=str(Peptquant_container.number_input(label="MAXIMUM:",min_value=-10,max_value=10,step=1,value=10))
        config_dict["LOGp_SNR_CUTOFF"]=str(Peptquant_container.number_input(label="LOG2_SNR_CUTOFF:",min_value=-5,max_value=5,step=1,value=1))
        protquant_expander=form.expander("Parameters for SIP protein quantification")
        config_dict["MIN_PEPTIDE_NUMBER"]=str(protquant_expander.number_input(label="MIN_PEPTIDE_NUMBER:",min_value=1,max_value=5,step=1,value=1))
        config_dict["MAX_CI_WIDTH"]=str(protquant_expander.number_input(label="MAX_CI_WIDTH:",min_value=1,max_value=10,step=1,value=5))
        config_dict["SMOOTHING_PROBABILITY_SPACE"]=str(protquant_expander.number_input(label="SMOOTHING_PROBABILITY_SPACE:",min_value=0.05,max_value=0.95,step=0.05,value=0.15))
        protquant_container=protquant_expander.container()
        protquant_container.write("")
        protquant_container.write("Log2 ratio for protein quantification:")
        config_dict["LOGT_MINIMUM"]=str(protquant_container.number_input(label="MINIMUM:",min_value=-10,max_value=10,step=1,value=-5))
        config_dict["LOGT_MAXIMUM"]=str(protquant_container.number_input(label="MAXIMUM:",min_value=-10,max_value=10,step=1,value=5))
        config_dict["MAX_LOG2_SNR"]=str(protquant_container.number_input(label="MAX_LOG2_SNR:",min_value=-5,max_value=5,step=1,value=4))
        config_dict["LOG2_RATIO_DISCRETIZATION"]=str(protquant_container.number_input(label="LOG2_RATIO_DISCRETIZATION:",min_value=0.0,max_value=1.0,step=0.1,value=0.1))

    config_dict["REMOVE_AMBIGUOUS_PEPTIDES"]  = quant_expander.selectbox(label='REMOVE_AMBIGUOUS_PEPTIDES: ',options=["true", "false"])

    submit_button = form.button(label='Run Process')
    return submit_button, config_dict

def drawCol2(col2:st,session_dict:dict,user_data_path:str,aval_job:str)->None:
    new_col = col2.container()
    tmp_label=functions.readUploadFile(os.path.join(user_data_path,session_dict['user_email']),"loading")
    if tmp_label=="0":
        df=pd.DataFrame(["No file"],columns=["Files"])
    else:
        result=set()
        for fname in os.listdir(os.path.join(user_data_path,session_dict['user_email'],aval_job)):
            if fname.find(".FT1")!=-1:
                result.add(fname.replace(".FT1",".raw"))
            elif fname.find(".FT2")!=-1:
                result.add(fname.replace(".FT2",".raw"))
            elif fname.find(".fasta") !=-1 or fname.find(".fa")!=-1:
                result.add(fname)
        df=pd.DataFrame(result,columns=["Files"])

    hide_table_row_index = """
            <style>      
            table {width: 100%;}
            thead tr th:first-child {display:none}
            tbody th {display:none}
            tr:nth-child(even){background-color: #f2f2f2; color:black;}
            tr:nth-child(odd){background-color: white; color:black;}
            th{
                background-color: #04AA6D;
            }
            </style>
            """
    new_col.markdown(hide_table_row_index, unsafe_allow_html=True)
    new_col.table(df)

def Sipros_Section():
    # session_dict={}
    st.markdown("""
<style>div[data-testid="stToolbar"] { display: none;}</style>
""", unsafe_allow_html=True)
    max_jobs = 1
    table_display_lable = True
    try:
        session_dict = st.session_state.session_infro 
    except:
        st.session_state.current_state='Landing Page'
        Landing_Page()
        return
    
    if session_dict['user_email'] not in os.listdir(user_data_path):
        os.makedirs(os.path.join(user_data_path,session_dict['user_email']))
        functions.Generateuserinfo(max_jobs,os.path.join(user_data_path,session_dict['user_email']),session_dict['user_email'])
        functions.generateUploadFile(os.path.join(user_data_path,session_dict['user_email']))
    session_dict['user_info'] = functions.Readuserinfo(os.path.join(user_data_path,session_dict['user_email']))

    aval_job = functions.Checkavab(session_dict['user_info'])
    if "session_infro" not in st.session_state:
        st.write("authorization expire!\n Please login again!\n")
        st.session_state.current_state='Landing Page'
    if aval_job != "-1" and  "0" in os.listdir(os.path.join(user_data_path,session_dict['user_email'])) and functions.readUploadFile(os.path.join(user_data_path,session_dict['user_email']),"loading")=="0":
        shutil.rmtree(os.path.join(user_data_path,session_dict['user_email'],"0"))
    if aval_job == "-1":
        st.write("The job queue is full. Waiting for the job completing!\nYou can close the browser.")
    elif aval_job != "-1" and table_display_lable: 
        config_dict = dict()  
        side_bar=st.sidebar
        side_bar.markdown("""<style>div[data-testid="stToolbar"] { display: none;} 
            div.css-6qob1r.e1fqkh3o3{ text-align:center;}
            div.stButton button:nth-child(1){border: none;height:100%; width:100%;}
            div.stButton button:nth-child(1):hover {color:#6495ED;}
            div.stButton button:nth-child(1):focus {border: none; color: #6495ED;box-shadow: 0.1px 0.11px 1px 1px #6495ED;}
            div.stButton button:nth-child(1):active {border: none; color: white;box-shadow: 0.1px 0.11px 1px 1px #6495ED;background-color:#6495ED}         
        </style>""", unsafe_allow_html=True)
        action1=side_bar.button("HOME")
        side_bar.write("")
        action2=side_bar.button("INSTRUCTION") 
        side_bar.write("")
        action3=side_bar.button("PUBLICATION")
        # side_bar.write("")
        # action6=side_bar.button("DEMO")  
        side_bar.write("")
        action4=side_bar.button("SOURCE CODE") 
        side_bar.write("")
        action5=side_bar.button("CONTACT US") 

        body_container=st.empty()
        container2=body_container.container()
        if action2:
            instructionPage(container2)
        if action3:
            publicationPage(container2)
        if action4:
            sourceCodePage(container2)
        if action5:
            contactPage(container2)
        # if action6:
        #     demoPage(action6, container2)
        if not action2 and not action3 and not action4 and not action5:
            col1,col2=container2.columns([3,1])
            col2=col2.empty()
        

            page_drive,page_params = col1.tabs(["Google Drive","Parameters"])
            page_drive.markdown("""
            <style> div[class="css-6kekos edgvbvh5"] {background-color:#6495ED;}
            </style>
            """, unsafe_allow_html=True)
            with page_drive.form(key="Floder link"):
                config_dict["folderlink"]= st.text_input('Enter google drive folder share link from the google drive folder (steps shown below) containing the following: raw files, and a protein database without decoy')
                drive_button=st.form_submit_button("Upload")
            page_drive.write("")
            page_drive.write("")
            page_drive.write("")
            page_drive.write("")
            page_drive.write("")
            page_drive.write("")
            page_drive.write("")
            page_drive.write("")
            page_drive.write("")
            page_drive.write("")
            page_drive.write("")
            page_drive.write("")
            page_drive.write("")
            page_drive.write("")
            page_drive.write("")
            page_drive.write("")
            page_drive.write("")
            page_drive.write("")
            page_drive.write("")
            page_drive.write("")
            page_drive.write("")

            drawCol2(col2,session_dict,user_data_path,aval_job)
            session_dict['whole_path'] = os.path.join(user_data_path,session_dict['user_email'],aval_job)
            if drive_button and config_dict["folderlink"]=="":
                page_drive.error("Upload your Google Drive Folder!\n")
            elif drive_button and config_dict["folderlink"]!="":
                logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.INFO)
                if "0" in os.listdir(os.path.join(user_data_path,session_dict['user_email'])):
                    shutil.rmtree(session_dict['whole_path'])
                os.mkdir(session_dict['whole_path'])
                st.session_state.page_drive = page_drive.empty()
                mul_file_label, err_label =Downloadfromdrive(session_dict,user_data_path, config_dict,page_drive.empty())
                functions.writeUploadFile(os.path.join(user_data_path,session_dict['user_email']),("multiple",str(int(mul_file_label))),"")
                if err_label == "miss file":
                    page_drive.error("You miss some raw/fasta files or system not support this raw file.\n Please upload the new raw and fasta files.")
                elif err_label == "no file":
                    page_drive.error("Your uploading folder is not in your google drive.\n Please upload the new folder.")
                else:
                    functions.writeUploadFile(os.path.join(user_data_path,session_dict['user_email']),("loading","1"),"")
                    drawCol2(col2,session_dict,user_data_path,aval_job)
            running_button=False
            running_button, config_dict=paramsTable(page_params,config_dict)
            tmp_label=functions.readUploadFile(os.path.join(user_data_path,session_dict['user_email']),"loading")
            if tmp_label=="1" and running_button:
                if True:
                    for fname in os.listdir(os.path.join(user_data_path,session_dict['user_email'],aval_job)):
                        if fname.find(".fa")!=-1:
                            functions.Reversedatabase(session_dict['whole_path'] +"/"+fname,session_dict['whole_path']+"/newdatabase.fasta",config_dict)
                            os.remove(session_dict['whole_path']+"/"+fname)
                    genconfigfiles.generateSiprosConfigFile(config_dict,session_dict["whole_path"])
                    genconfigfiles.generateProrataConfigFile(config_dict,session_dict["whole_path"])
                    if config_dict["Version_selection"]=="single machine":
                        functions.Generateconfigfile(session_dict,3.0,False,config_dict["Search_Type"])
                    else:
                        functions.Generateconfigfile(session_dict,3.0,True,config_dict["Search_Type"])

                    functions.Runonserver(host_inf,session_dict,False,"sbatchscript.sh", True)
                    logging.info("job id:"+str(session_dict['job_id']))
                    session_dict['user_info'][session_dict['whole_path'].split("/")[-1]][0]="1"
                    functions.Writeuserinfo(session_dict,os.path.join(user_data_path,session_dict['user_email']))
                    body_container.write("Running!\n The results will be sent back to your Gmail. You can close the browser.")
                    time.sleep(5)
                    statues_running = ""
                    while True:
                        try:
                            check_status, running_time = functions.Checkstatus(host_inf,session_dict['job_id'],session_dict["whole_path"])
                        except:
                            print("sever error!")
                            check_status="RUNNING"
                            running_time=""
                    
                        if check_status == "COMPLETED":
                            statues_running="COMPLETED"
                            break
                        elif check_status == "PENDING" or check_status == "RUNNING" or check_status == "REQUEUED" or check_status == "RESIZING" or check_status == "REVOKED" or check_status == "SUSPENDED":
                            logging.info('watting')
                        elif check_status == "CANCELLED+" or check_status == "NODE_FAIL" or check_status == "PREEMPTED": 
                            statues_running = "FAILED"
                            logging.info('CANCELLED')
                            break
                        elif check_status == "DEADLINE" or check_status == "TIMEOUT":
                            statues_running = "deadline"
                            logging.info('DEADLINE')
                            if running_time.find("-")!=-1:
                                if running_time=="1-23:59:59" or running_time.split("-")[0]=="2":#hours
                                    break
                            if config_dict["Version_selection"]=="single machine":
                                functions.Generateconfigfile(session_dict,3.0,False,config_dict["Search_Type"])
                            else:
                                if not functions.Generateconfigfile(session_dict,3.0,True,config_dict["Search_Type"]):
                                    config_dict["Version_selection"]="single machine"
                            functions.Runonserver(host_inf,session_dict,True,"sbatchscript.sh",True)
                        elif check_status == "BOOT_FAIL":
                            statues_running = "FAILED"
                            logging.info('BOOT_FAIL')
                            break
                        elif check_status == "FAILED":
                            statues_running = "FAILED"
                            logging.info('FAILED')
                            break
                        elif check_status == "OUT_OF_MEMORY":
                            statues_running = "OUT_OF_MEMORY"
                            logging.info('OUT_OF_MEMORY')
                            break
                        else:
                            statues_running = "FAILED"
                            logging.info('else')
                            break
                        st.session_state.session_infro = session_dict
                        time.sleep(600)
                    if config_dict["Version_selection"]!="single machine" and statues_running == "COMPLETED":
                        functions.copyFiles(host_inf,session_dict,config_dict["Search_Type"])
                        functions.Runonserver(host_inf,session_dict,False,"sbatchscript2.sh", False)
                        while True:
                            try:
                                check_status, running_time = functions.Checkstatus(host_inf,session_dict['job_id'],session_dict["whole_path"])
                            except:
                                print("server error!")
                                check_status="RUNNING"
                                running_time=""
                            if check_status == "COMPLETED":
                                statues_running="COMPLETED"
                                break
                            elif check_status == "PENDING" or check_status == "RUNNING" or check_status == "REQUEUED" or check_status == "RESIZING" or check_status == "REVOKED" or check_status == "SUSPENDED":
                                logging.info('watting')
                            elif check_status == "CANCELLED+" or check_status == "NODE_FAIL" or check_status == "PREEMPTED": 
                                statues_running = "FAILED"
                                logging.info('CANCELLED')
                                break
                            elif check_status == "DEADLINE" or check_status == "TIMEOUT":
                                statues_running = "deadline"
                                logging.info('DEADLINE')
                                if running_time.find("-")!=-1:
                                    if running_time=="1-23:59:59" or running_time.split("-")[0]=="2":#hours
                                        break
                                if config_dict["Version_selection"]=="single machine":
                                    functions.Generateconfigfile(session_dict,3.0,False,config_dict["Search_Type"])
                                else:
                                    if not functions.Generateconfigfile(session_dict,3.0,True,config_dict["Search_Type"]):
                                        config_dict["Version_selection"]="single machine"
                                functions.Runonserver(host_inf,session_dict,False,"sbatchscript2.sh", False) 
                            elif check_status == "BOOT_FAIL":
                                statues_running = "FAILED"
                                logging.info('BOOT_FAIL')
                                break
                            elif check_status == "FAILED":
                                statues_running = "FAILED"
                                logging.info('FAILED')
                                break
                            elif check_status == "OUT_OF_MEMORY":
                                statues_running = "OUT_OF_MEMORY"
                                logging.info('OUT_OF_MEMORY')
                                break
                            else:
                                statues_running = "FAILED"
                                logging.info('FAILED')
                                break
                            st.session_state.session_infro = session_dict
                            time.sleep(600)
                    if statues_running == "COMPLETED":
                        functions.Downloadfile(host_inf,session_dict,config_dict["Search_Type"])           
                        email_massge = " Dear Sipros user,\n Your job completed.\n"
                        try:
                            functions.sendEmail(email_massge,session_dict,True)
                        except Exception as e:
                            print(e)
                        else:
                            functions.delOnSlurm(host_inf,session_dict)
                        session_dict['user_info'][session_dict['whole_path'].split("/")[-1]][0]="0"
                        functions.Writeuserinfo(session_dict,os.path.join(user_data_path,session_dict['user_email']))    
                        functions.generateUploadFile(os.path.join(user_data_path,session_dict['user_email']))
                        col1.success("Finish!")           
                    elif statues_running == "FAILED":
                        col1.error("FAILED")
                        email_massge = "Dear Sipros user,\n Your job not completed. Please check your Sipros parameters, raw files, and database.\n"
                        try:
                            functions.sendEmail(email_massge,session_dict,False)
                        except Exception as e:
                            print(e)
                        else:
                            functions.delOnSlurm(host_inf,session_dict)
                        session_dict['user_info'][session_dict['whole_path'].split("/")[-1]][0]="0"
                        functions.Writeuserinfo(session_dict,os.path.join(user_data_path,session_dict['user_email']))
                        functions.generateUploadFile(os.path.join(user_data_path,session_dict['user_email']))
                    elif statues_running == "deadline":
                        col1.error("running time out of limitation!")
                        email_massge = "Dear Sipros user,\n Your job is out of time limitation. Please submit samller job.\n"
                        try:
                            functions.sendEmail(email_massge,session_dict,False)
                        except Exception as e:
                            print(e)
                        else:
                            functions.delOnSlurm(host_inf,session_dict)
                        session_dict['user_info'][session_dict['whole_path'].split("/")[-1]][0]="0"
                        functions.Writeuserinfo(session_dict,os.path.join(user_data_path,session_dict['user_email']))
                        functions.generateUploadFile(os.path.join(user_data_path,session_dict['user_email']))
                    elif statues_running =="OUT_OF_MEMORY":
                        col1.error("memory out of limitation!")
                        email_massge = "Dear Sipros user,\n Your job is out of memory. Please submit samller job.\n"
                        try:
                            functions.sendEmail(email_massge,session_dict,False)
                        except Exception as e:
                            print(e)
                        else:
                            functions.delOnSlurm(host_inf,session_dict)
                        session_dict['user_info'][session_dict['whole_path'].split("/")[-1]][0]="0"
                        functions.Writeuserinfo(session_dict,os.path.join(user_data_path,session_dict['user_email']))
                        functions.generateUploadFile(os.path.join(user_data_path,session_dict['user_email']))
                    # Sign_out()
        bottom_container=st.container()
        bottom_container.markdown('''
        <div style="position:relative;right:0px;bottom:0px;height:200px;background-color:#6495ED;">
            <div style="text-align:center; position:relative; top:4px">
                <p> Developed by OU and UNT</p>
            </div>
            <div style="margin: auto; width:110px;">
                <img src="data:image/png;base64,{0}" alt="OU" style="border-radius:10%;  height:50px; width:50px;">
                <img src="data:image/png;base64,{1}" alt="UNT" style="border-radius:10%;   height:50px; width:50px;">
            </div>
            <div style="text-align:center; position:relative; top:4px">
                <p> Funded by NCCIH and NIGMS.</p>
            </div>
            <div style="margin: auto; width:110px;">
                <img src="data:image/png;base64,{2}" alt="OU" style="border-radius:10%;  height:50px; width:50px;">
                <img src="data:image/png;base64,{3}" alt="UNT" style="border-radius:10%;   height:50px; width:50px;">
            </div>
        </div>'''.format(base64.b64encode(Path("resources/University_of_Oklahoma_Logo.png").read_bytes()).decode(),base64.b64encode(Path("resources/UNT_Logo_Alternative2.png").read_bytes()).decode(),base64.b64encode(Path("resources/nccih.png").read_bytes()).decode(),base64.b64encode(Path("resources/nigms.jpeg").read_bytes()).decode()),unsafe_allow_html=True)
    
    return


def Landing_Page():
    st.markdown("""<style>div[data-testid="stToolbar"] { display: none;} 
            div.css-6qob1r.e1fqkh3o3{ text-align:center;}
            div.stButton>button:nth-child(1){border: none;height:100%; width:100%;}
            div.stButton>button:nth-child(1):hover {color:#6495ED;}
            div.stButton>button:nth-child(1):focus {border: none; color: #6495ED;box-shadow: 0.1px 0.11px 1px 1px #6495ED;}
            div.stButton>button:nth-child(1):active {border: none; color: white;box-shadow: 0.1px 0.11px 1px 1px #6495ED;background-color:#6495ED}         
        </style>""", unsafe_allow_html=True)
    side_bar = st.sidebar
    container=st.container()

    action1=side_bar.button("HOME")
    side_bar.write("")
    action2=side_bar.button("INSTRUCTION") 
    side_bar.write("")
    action3=side_bar.button("PUBLICATION") 
    # side_bar.write("")
    # action6=side_bar.button("DEMO") 
    side_bar.write("")
    action4=side_bar.button("SOURCE CODE") 
    side_bar.write("")
    action5=side_bar.button("CONTACT US")



    body_container=container.container()
    col21,col22,col23=body_container.columns([1,5,1])

    col22_empty=col22.empty()
    login_container=col22_empty.container()
    if not action1 and not action2 and not action3 and not action4 and not action5:
        loginPage(login_container)
    if action1:
        loginPage(login_container)
    if action2:
        instructionPage(login_container)
    if action3:
        publicationPage(login_container)
    if action4:
        sourceCodePage(login_container)
    if action5:
        contactPage(login_container)
    # if action6:
    #     paramPage(login_container)
    bottom_container=container.container()
    bottom_container.markdown('''
        <div style="position:relative;right:0px;bottom:0px;height:200px;background-color:#6495ED;">
            <div style="text-align:center; position:relative; top:4px">
                <p> Developed by OU and UNT</p>
            </div>
            <div style="margin: auto; width:110px;">
                <img src="data:image/png;base64,{0}" alt="OU" style="border-radius:10%;  height:50px; width:50px;">
                <img src="data:image/png;base64,{1}" alt="UNT" style="border-radius:10%;   height:50px; width:50px;">
            </div>
            <div style="text-align:center; position:relative; top:4px">
                <p> Funded by NCCIH, NIGMS, and NLM.</p>
            </div>
            <div style="margin: auto; width:160px;">
                <img src="data:image/png;base64,{2}" alt="OU" style="border-radius:10%;  height:50px; width:50px;">
                <img src="data:image/png;base64,{3}" alt="UNT" style="border-radius:10%;   height:50px; width:50px;">
                <img src="data:image/png;base64,{4}" alt="UNT" style="border-radius:10%;   height:50px; width:50px;">
            </div>
        </div>'''.format(base64.b64encode(Path("resources/University_of_Oklahoma_Logo.png").read_bytes()).decode(),base64.b64encode(Path("resources/UNT_Logo_Alternative2.png").read_bytes()).decode(),base64.b64encode(Path("resources/nccih.png").read_bytes()).decode(),base64.b64encode(Path("resources/nigms.jpeg").read_bytes()).decode(),base64.b64encode(Path("resources/nlm.png").read_bytes()).decode()),unsafe_allow_html=True)
    # if action6:
    #     demoPage(action6,login_container)

def instructionPage(col22:st)->None:
    col22.write("")
    col22.write("")
    col22.header("Instruction:")
    col22.subheader("The steps:")
    col22.write("1. Upload your raw data files and database file into one folder on Google Drive.")
    img1=Image.open("resources/folder.png")
    col22.image(img1,width=400)
    col22.write("2. Get link of data folder.")
    col22.caption(" 2.1 Right click the folder and select get link.")
    img2=Image.open("resources/get_link.png")
    col22.image(img2,width=400)
    col22.caption(" 2.2 Click copy link button.")
    img3=Image.open("resources/copy_link.png")
    col22.image(img3,width=400)
    col22.write("3. Login with Google account.")
    col22.caption(" 3.1 Click Login button.")
    img4=Image.open("resources/login.png")
    col22.image(img4,width=400)
    col22.caption(" 3.2 Select Google account.")
    img5=Image.open("resources/select_account.png")
    col22.image(img5,width=400)
    col22.caption(" 3.3 Authorize Gmail and Google Drive.")
    img6=Image.open("resources/Picture1.png")
    col22.image(img6,width=400)
    img7=Image.open("resources/Picture2.png")
    col22.image(img7,width=400)
    img8=Image.open("resources/authorize.png")
    col22.image(img8,width=400)
    col22.write("4. Upload data to server.")
    col22.caption(" Paste copied link and click upload button. Waiting for finishing to upload files. ")
    img9=Image.open("resources/google_link.png")
    col22.image(img9,width=400)
    col22.write("5. Set paramerters.")
    col22.caption(" After program run, close the browser. The results will be sent back by eamil.")
    img10=Image.open("resources/parameters.png")
    col22.image(img10,width=400)
    col22.subheader("Some Parameters:")
    col22.write("")
    df=pd.DataFrame(ptm_expanation)

    hide_table_row_index = """
            <style>      
            table {width: 100%;}
            thead tr th:first-child {display:none}
            tbody th {display:none}
            tr:nth-child(even){background-color: #f2f2f2; color:black;}
            tr:nth-child(odd){background-color: white; color:black;}
            th{
                background-color: #04AA6D;
            }
            </style>
            """
    col22.markdown(hide_table_row_index, unsafe_allow_html=True)
    col22.table(df)







def contactPage(col22:st)->None:
    col22.markdown("""
    <style>div[class="css-6vw1hb e1tzin5v2"] {height: 700px;}
    .block-container div[class="css-1mapwt e1tzin5v0"] div[class="css-1mapwt e1tzin5v0"] {height: 900px;}
    </style>
    """,unsafe_allow_html=True)
    col22.write("")
    col22.write("")
    col22.header("Contact Imformation:")
    col22.write("")
    col22.write("""<p>Chongle Pan</p><p><div style="display: inline-block;"><b>Email:</b></div><div style="display: inline-block;">cpan@ou.edu</div></p>""",unsafe_allow_html=True)
    col22.write("")
    col22.write("""<p>Xuan Guo</p><p><div style="display: inline-block;"><b>Email:</b></div><div style="display: inline-block;">xuan.guo@unt.edu</div></p>""",unsafe_allow_html=True)
    col22.write("")
    col22.write("""<p>Jiancheng Li</p><p><div style="display: inline-block;"><b>Email:</b></div><div style="display: inline-block;">jianchengli@my.unt.edu</div></p>""",unsafe_allow_html=True)



def publicationPage(col22:st)->None:
    url1="https://pubs.acs.org/doi/full/10.1021/ac060654b"
    url2="https://academic.oup.com/bioinformatics/article/34/5/795/4209993"
    url3="https://academic.oup.com/bioinformatics/article/29/16/2064/200485"
    url4="https://academic.oup.com/bioinformaticsadvances/advance-article/doi/10.1093/bioadv/vbae024/7613684"
    col22.markdown("""
    <style>div[class="css-6vw1hb e1tzin5v2"] {height: 700px;} 
    div[class="css-1mapwt e1tzin5v0"] {height: 900px;}
    </style>
    """,unsafe_allow_html=True)
    col22.write("")
    col22.write("")
    col22.header("Publications:") 
    col22.write(f'''<p>1. Li, J., Xiong, Y., Feng, S., Pan, C., & Guo, X. (2024). CloudProteoAnalyzer: scalable processing of big data from proteomics using cloud computing.<a href={url4}> Bioinformatics Advances, vbae024</a></p>''',unsafe_allow_html=True)
    col22.write(f'''<p>2. Guo, X., Li, Z., Yao, Q., Mueller, R.S., Eng, J.K., Tabb, D.L., Hervey IV, W.J. and Pan, C., 2018. Sipros ensemble improves database searching and filtering for complex metaproteomics.<a href={url2}> Bioinformatics, 34(5), pp.795-802</a></p>''',unsafe_allow_html=True)
    col22.write(f'''<p>3. Wang, Y., Ahn, T.H., Li, Z. and Pan, C., 2013. Sipros/ProRata: a versatile informatics system for quantitative community proteomics.<a href={url3}> Bioinformatics, 29(16), pp.2064-2065</a></p>''',unsafe_allow_html=True)
    col22.write(f'''<p>4. Pan, C., Kora, G., McDonald, W.H., Tabb, D.L., VerBerkmoes, N.C., Hurst, G.B., Pelletier, D.A., Samatova, N.F. and Hettich, R.L., 2006. ProRata: a quantitative proteomics program for accurate protein abundance ratio estimation with confidence interval evaluation.<a href={url1}> Analytical chemistry, 78(20), pp.7121-7131</a></p''',unsafe_allow_html=True)


def sourceCodePage(col22:st)->None:
    # col22=col22.container()
    col22.markdown("""
    <style>div[class="css-6vw1hb e1tzin5v2"] {height: 700px;} 
        div[class="css-cg3mbs e1tzin5v0"]  div[class="css-cg3mbs e1tzin5v0"] {height: 700px;}
    </style>
    """,unsafe_allow_html=True)
    col22.write("")
    col22.write("")
    col22.subheader("You can find the source code from our github.")
    col22.write("The link:")
    col22.write("https://github.com/Biocomputing-Research-Group/CloudProteoAnalyzer")
    col22.write("")
    col22.write("")
    col22.write("")
    col22.write("")
    col22.write("")
    col22.write("")
    col22.write("")
    col22.write("")
    col22.write("")
    col22.write("")
    col22.write("")
    col22.write("")
    col22.write("")
    col22.write("")
    col22.write("")
    col22.write("")
    return


def paramPage(col22:st)->None:
    col22.markdown("""
    <style>div[class="css-6vw1hb e1tzin5v2"] {height: 700px;} 
        div[class="css-cg3mbs e1tzin5v0"]  div[class="css-cg3mbs e1tzin5v0"] {height: 700px;}
    </style>
    """,unsafe_allow_html=True)
    col22.write("")
    col22.write("")
    col22.title("Guide to Sipros/Prorata results")
    col22.write("1.$/{output_dir/}$/{ms2_filename/}_$/{search_name/}.psm.txt : This file contains PSM level results from the PSM filtering. The columns in this file are")


def loginPage(col22:st)->None:
    flow = google_auth_oauthlib.flow.Flow.from_client_secrets_file(
        google_auth_file,
        scopes = ['https://www.googleapis.com/auth/drive','https://www.googleapis.com/auth/userinfo.email'])  
    flow.redirect_uri = web_url
    authorization_url, state = flow.authorization_url(
        # Enable offline access so that you can refresh an access token without
        # re-prompting the user for permission. Recommended for web server apps.
        access_type='offline',
        # Enable incremental authorization. Recommended as a best practice.
        include_granted_scopes='true',
        # prompt='consent'
    )
    url = authorization_url
        # url = 'https://stackoverflow.com'
    col22.markdown("""
    <style>div[class="css-6vw1hb e1tzin5v2"] {text-align:center;} </style>
    """,unsafe_allow_html=True)
    col22.title('Welcome to CloudProteoAnalyzer')
    col22.write('')
    col22.write('')
    col22.write('') 
    col22.write('Sign in with Google with the button below')
    col22.write(f'''<p><a href={url}><button style="background-color:#6495ED;border-radius: 8px;color:white;">Log in</button></a></p>''',unsafe_allow_html=True)
    data_path="/home/jiancheng/example"
    col22.write("")
    col22.write("")
    col22.write("")
    col22.write("")
    col22.write("")
    col22.write("")
    col22.write("")
    col22.title("DEMO PART")
    condition=col22.button("Running Demo Data")
    hide_table_row_index = """
            <style>      
            table {width: 100%;}
            thead tr th:first-child {display:none}
            tbody th {display:none}
            tr:nth-child(even){background-color: #f2f2f2; color:black;}
            tr:nth-child(odd){background-color: white; color:black;}
            th{
                background-color: #04AA6D;
            }
            </style>
            """
    if not condition:
        result=set()
        for fname in os.listdir(data_path):
            if fname.find(".FT1")!=-1:
                result.add(fname.replace(".FT1",".raw"))
            elif fname.find(".fasta") !=-1 or fname.find(".fa")!=-1:
                result.add(fname)
        df=pd.DataFrame(result,columns=["Files"])  
        col22.markdown(hide_table_row_index, unsafe_allow_html=True)
        col22.table(df)         
    else:
        result=set()
        for fname in os.listdir(data_path):
            if fname.find(".FT1")!=-1:
                result.add(fname.replace(".FT1",".raw"))
            elif fname.find(".fasta") !=-1 or fname.find(".fa")!=-1:
                result.add(fname)
            elif fname.find(".txt")!=-1:
                result.add(fname)
            elif fname.find(".tab")!=-1:
                result.add(fname)
        df=pd.DataFrame(result,columns=["Files"])
        col22.markdown(hide_table_row_index, unsafe_allow_html=True)
        col22.table(df)
    expan=col22.expander("Guide to Sipros/Prorata results")
    expan.markdown(''' <style>div[data-testid="stExpander"] {text-align:left;} </style>''', unsafe_allow_html=True)

    expan.subheader("1.{ms2_filename}_{search_name}.psm.txt : This file contains PSM level results from the PSM filtering. The columns in this file are:")
    expan.write("   Filename = Filename of input MS2 file")
    expan.write("   ScanNumber = Scan number of the PSM")
    expan.write("   ParentCharge = Charge state of the PSM")
    expan.write("   MeasuredParentMass = Measured parent mass")
    expan.write("   CalculatedParentMass = Calculated parent mass from peptide sequence")
    expan.write("   MassErrorDa = Mass error in Da with 1-Da error correction")
    expan.write("   MassErrorPPM = Mass error in PPM with 1-Da error correction")
    expan.write("   ScanType = Scan type of the PSM")
    expan.write("   SearchName = Sipros search name")
    expan.write("   ScoringFunction = Scoring function used in the search")
    expan.write("   Score = Predicted Probability of being true PSM")
    expan.write("   DeltaZ = Difference between the best PSM score and the next best PSM of this scan")
    expan.write("   DeltaP = Difference between the best modified PSM and its PTM isoform")
    expan.write("   IdentifiedPeptide = Identified peptide sequence with potential PTMs and mutations")
    expan.write("   OriginalPeptide = Original peptide sequence in the FASTA file")
    expan.write("   ProteinNames = Names of proteins of the peptide")
    expan.write("   ProteinCount = Number of proteins that the peptide can be assigned to")
    expan.write("   TargetMatch = T for target match and F for decoy match")
    expan.subheader("2.{ms2_filename}_{search_name}.pep.txt : This file contains Peptide level results from the Peptide filtering. The columns in this file are:")
    expan.write("   IdentifiedPeptide = Identified peptide sequence with potential PTMs and mutations")
    expan.write("   ParentCharge = Charge state of identified peptide")
    expan.write("   OriginalPeptide = Original peptide sequence in the FASTA file")
    expan.write("   ProteinNames = Names of proteins of the peptide")
    expan.write("   ProteinCount = Number of proteins that the peptide can be assigned to")
    expan.write("   TargetMatch = T for target match and F for decoy match")
    expan.write("   SpectralCount = Number of PSMs in which the peptide is identified")
    expan.write("   BestScore = The best score of those PSMs")
    expan.write("   PSMs = List of PSMs for the peptide: MS2_Filename[Scan_Number]")
    expan.write("   ScanType = Scan type of those PSMs")
    expan.write("   SearchName = Sipros search name")
    expan.subheader("3.{ms2_filename}_{search_name}.pro.txt : This file contains Protein level results from the Protein assembling. The columns in this file are:")
    expan.write("   ProteinID = Names of the protein")
    expan.write("    Run#_UniquePeptideCounts = Number of unique peptides in a run")
    expan.write("   Run#_TotalPeptideCounts = Number of all peptides in a run")
    expan.write("   Run#_UniqueSpectrumCounts = Number of unique PSM in a run")
    expan.write("   Run#_TotalSpectrumCounts = Number of all PSM in a run")
    expan.write("   Run#_BalancedSpectrumCounts = Balanced spectrum count in a run")
    expan.write("   Run#_NormalizedBalancedSpectrumCounts = Normalized Balanced spectrum count in a run")
    expan.write("   ProteinDescription = Protein description")
    expan.write("   TargetMatch = T for target match and F for decoy match")
    expan.subheader("4.{ms2_filename}_{search_name}.pro2psm.txt : This file contains the spectrum count and related statistics for each identified protein. The columns in this file are:")
    expan.write("   + = Marker of a protein line")
    expan.write("   ProteinID = Names of the protein")
    expan.write("   Run#_UniquePeptideCounts = Number of unique peptides in a run")
    expan.write("   Run#_TotalPeptideCounts = Number of all peptides in a run")
    expan.write("   Run#_UniqueSpectrumCounts = Number of unique PSM in a run")
    expan.write("   Run#_TotalSpectrumCounts = Number of all PSM in a run")
    expan.write("   Run#_BalancedSpectrumCounts = Balanced spectrum count in a run")
    expan.write("   Run#_NormalizedBalancedSpectrumCounts = Normalized Balanced spectrum count in a run")
    expan.write("   ProteinDescription = Protein description")
    expan.write("   TargetMatch = T for target match and F for decoy match")
    expan.write("   * = Marker of a peptide line")
    expan.write("   IdentifiedPeptide = Identified peptide sequence with potential PTMs and mutations")
    expan.write("   ParentCharge = Charge state of identified peptide")
    expan.write("   OriginalPeptide = Original peptide sequence in the FASTA file")
    expan.write("   ProteinNames = Names of proteins of the peptide")
    expan.write("   ProteinCount = Number of proteins that the peptide can be assigned to")
    expan.write("   TargetMatch = T for target match and F for decoy match")
    expan.write("   SpectralCount = Number of PSMs in which the peptide is identified")
    expan.write("   BestScore = The best score of those PSMs")
    expan.write("   PSMs = List of PSMs for the peptide: MS2_Filename[Scan_Number]")
    expan.write("   ScanType = Scan type of those PSMs")
    expan.write("   SearchName = Sipros search name")
    expan.subheader("5.{ms2_filename}_{search_name}.pro2pep.txt : This file contains the peptide count and related statistics for each identified protein. The columns in this file are")
    expan.write("   + = Marker of a protein line")
    expan.write("   ProteinID = Names of the protein")
    expan.write("   Run#_UniquePeptideCounts = Number of unique peptides in a run")
    expan.write("   Run#_TotalPeptideCounts = Number of all peptides in a run")
    expan.write("   Run#_UniqueSpectrumCounts = Number of unique PSM in a run")
    expan.write("   Run#_TotalSpectrumCounts = Number of all PSM in a run")
    expan.write("   Run#_BalancedSpectrumCounts = Balanced spectrum count in a run")
    expan.write("   ProteinDescription = Protein description")
    expan.write("   TargetMatch = T for target match and F for decoy match")
    expan.write("   * = Marker of a PSM line")
    expan.write("   Filename = Filename of input MS2 file")
    expan.write("   ScanNumber = Scan number of the PSM")
    expan.write("   ParentCharge = Charge state of the PSM")
    expan.write("   MeasuredParentMass = Measured parent mass")
    expan.write("   CalculatedParentMass = Calculated parent mass from peptide sequence")
    expan.write("   MassErrorDa = Mass error in Da with 1-Da error correction")
    expan.write("   MassErrorPPM = Mass error in PPM with 1-Da error correction")
    expan.write("   ScanType = Scan type of the PSM")
    expan.write("   SearchName = Sipros search name")
    expan.write("   ScoringFunction = Scoring function used in the search")
    expan.write("   Score = Score")
    expan.write("   DeltaZ = Difference between the best PSM and the next best PSM of this scan")
    expan.write("   DeltaP = Difference between the best modified PSM and its PTM isoform")
    expan.write("   IdentifiedPeptide = Identified peptide sequence with potential PTMs and mutations")
    expan.write("   OriginalPeptide = Original peptide sequence in the FASTA file")
    expan.write("   ProteinNames = Names of proteins of the peptide")
    expan.write("   ProteinCount = Number of proteins that the peptide can be assigned to")
    expan.write("   TargetMatch = T for target match and F for decoy match")
    expan.subheader("6.{ms2_filename}_{search_name}_Spe2Pep.txt: This file is an intermediate file generated by the database-searching of Sipros Ensemble. The columns in this file are:")
    expan.write("   + = Marker of a spectrum line")
    expan.write("   Filename = Filename of input MS2 file")
    expan.write("   ScanNumber = Scan number of the PSM")
    expan.write("   ParentCharge = Charge state of the PSM")
    expan.write("   MeasuredParentMass = Measured parent mass")
    expan.write("   ScanType = Scan type of the PSM")
    expan.write("   SearchName = Sipros search name")
    expan.write("   TotalIntensity = Sum of all the peak intensities of the PSM")
    expan.write("   MaxIntensity = Maximum intensity of the PSM")
    expan.write("   * = Marker of a peptide line")
    expan.write("   IdentifiedPeptide = Identified peptide sequence with potential PTMs and mutations")
    expan.write("   OriginalPeptide = Original peptide sequence in the FASTA file")
    expan.write("   CalculatedParentMass = Calculated parent mass from peptide sequence")
    expan.write("   MVH = Multivariate hypergeometric Score")
    expan.write("   Xcorr = Cross correlation score")
    expan.write("   WDP = Weighted dot product score")
    expan.write("   ProteinNames = Names of proteins of the peptide")
    expan.subheader("7.{ms2_filename}_{search_name}.tab : This file is an intermediate file generated by summarizing all the Spe2Pep.txt files. The columns are the features could be used for ensemble learning steps. The columns in this file are:")
    expan.write("   FileName = Filename of input MS2 file")
    expan.write("   ScanNumber = Scan number of the PSM")
    expan.write("   ParentCharge = Charge state of the PSM")
    expan.write("   MeasuredParentMass = Measured parent mass")
    expan.write("   ScanType = Scan type of the PSM")
    expan.write("   SearchName = Sipros search name")
    expan.write("   IdentifiedPeptide = Identified peptide sequence with potential PTMs and mutations")
    expan.write("   OriginalPeptide = Original peptide sequence in the FASTA file")
    expan.write("   CalculatedParentMass = Calculated parent mass from peptide sequence")
    expan.write("   MVH = Score by using the MVH scoreing function")
    expan.write("   Xcorr = Score by using the Xcorr scoring function") 
    expan.write("   WDP = Score by using the weighted-dot-product scoring function")
    expan.write("   ProteinNames = Names of proteins of the peptide")
    expan.write("   ScoreAgreement = Count of scores that rank the current PSM as the top one")
    expan.write("   DeltaRP1 = Fractional difference between current and best Rank Product based on MVH")
    expan.write("   DeltaRP2 = Fractional difference between current and best Rank Product based on XCorr")
    expan.write("   DeltaRP3 = Fractional difference between current and best Rank Product based on WDP")
    expan.write("   DeltaRS1 = Fractional difference between current and best MVH") 
    expan.write("   DeltaRS2 = Fractional difference between current and best XCorr")
    expan.write("   DeltaRS3 = Fractional difference between current and best WDP")
    expan.write("   DiffRP1 = Difference between current and next best Rank Product based on MVH")
    expan.write("   DiffRP2 = Difference between current and next best Rank Product based on Xcorr")
    expan.write("   DiffRP3 = Difference between current and next best Rank Product based on WDP")
    expan.write("   DiffRS1 = Difference between current and next best MVH")
    expan.write("   DiffRS2 = Difference between current and next best Xcorr")
    expan.write("   DiffRS3 = Difference between current and next best WDP")
    expan.write("   DiffNorRP1 = Fractional difference between current and next best Rank Product based on MVH ")
    expan.write("   DiffNorRP2 = Fractional difference between current and next best Rank Product based on Xcorr")
    expan.write("   DiffNorRP3 = Fractional difference between current and next best Rank Product based on WDP")
    expan.write("   DiffNorRS1 = Fractional difference between current and next best MVH")
    expan.write("   DiffNorRS2 = Fractional difference between current and next best Xcorr")
    expan.write("   DiffNorRS3 = Fractional difference between current and next best WDP")
    expan.write("   RetentionTime =  Retention time")
    expan.write("   LocalRank = Rank by using the rank product")
    expan.write("   DeltaP = Difference between the best modified PSM and its PTM isoform")
    col22.write("")
    col22.write("")
    col22.write("")
    col22.write("")


         

def Google_API_Section():
    g_session_dict= dict()

    os.environ['OAUTHLIB_RELAX_TOKEN_SCOPE'] = '1'
    if 'first_time_in_google_section' not in st.session_state:
        flow = google_auth_oauthlib.flow.Flow.from_client_secrets_file(google_auth_file,
            scopes = ['https://www.googleapis.com/auth/drive','https://www.googleapis.com/auth/userinfo.email'])

        flow.redirect_uri = web_url
        app_state = st.experimental_get_query_params()
        auth_code = app_state['code'][0]
        try:
            auth_token = flow.fetch_token(code=auth_code)
            user_info_service = gdiscovery.build(serviceName='oauth2', version='v2', credentials=flow.credentials)
            user_info = user_info_service.userinfo().get().execute()
            email_address = user_info["email"]
        except:
            # st.write("Google authentication expired!\n Please login again.")
            raise
            return
        else:
            g_session_dict['user_email']=email_address.split("@")[0]
            st.session_state.first_time_in_google_section = False
            st.session_state.oauth_flow = flow
            return g_session_dict
   


def Sign_out():
    #This section is responsible for signing out the user and destroying the oauth token
    #to revoke a token, a request needs to be made to "my_domain"/oauth/revoke    
    flow = st.session_state.oauth_flow
    credentials = flow.credentials
    sign_out_response = requests.post('https://oauth2.googleapis.com/revoke',
        params={'token': credentials.token},
        headers = {'content-type': 'application/x-www-form-urlencoded'})
    print("sign out")
    if sign_out_response.status_code == 200:
        st.success('Sign out complete. Thank you for using Sipros!')
    else:
        st.error('Unexpected response')
            #print("Unexpected response after attempted sign out. Response: ", sign_out_response)
    st.stop()
        
        

# +
 
state_dict = {
      'Landing Page': 0,
    #   'Table Section': 1,
      'Google Section': 1,
      'Sipros Section': 2,
      'Sign Out': 3
    }
app_state = st.experimental_get_query_params()
#On first iteration of initial window
if 'code' not in app_state and 'current_state' not in st.session_state:
    st.session_state.current_state = 'Landing Page'
    st.session_state.previous_state = 'Landing Page'
    st.session_state.google_section = True

#On first iteration of second window
elif 'code' in app_state and 'isFirst' not in st.session_state:
    try:
        st.session_state.session_infro = Google_API_Section()
    except:
        # st.write("Google authentication expired!\n Please login again.")
        st.session_state.current_state='Landing Page'
    st.session_state.current_state = 'Sipros Section'
    st.session_state.previous_state = 'Sipros Section'
    st.session_state.isFirst = False
    st.session_state.google_section = True
#     st.session_state.current_state = st.sidebar.selectbox('Statebox',['Landing Page',
#                                                   'Google Section','Sipros Section','Sign Out'], index=state_dict['Sipros Section'])#'Table Section','Sipros Section'

# else:
#     st.session_state.current_state = st.sidebar.selectbox('Statebox',['Landing Page',
                                                #   'Google Section','Sipros Section','Sign Out'],index=state_dict[st.session_state.previous_state])
if st.session_state.current_state == 'Landing Page':
    Landing_Page()
# elif st.session_state.current_state == 'Table Section':
#     #print("Entering Table Section")
#     Table_Section()
# elif st.session_state.current_state == 'Google Section':
#     #print("Entering Google Section")
#     Google_API_Section()
elif st.session_state.current_state == 'Sipros Section':
    #print("Entering Sipros Section")
    Sipros_Section()

#     try:
#         Sipros_Section()
#     except:
#         st.session_state.current_state='Landing Page'
elif st.session_state.current_state == 'Sign Out':
    #print("Entering Sign Out Section")
    Sign_out()
    
st.session_state.current_state = st.session_state.previous_state

# -


