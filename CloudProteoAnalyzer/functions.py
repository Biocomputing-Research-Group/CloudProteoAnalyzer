from distutils.command.config import config
from xmlrpc.client import boolean
import streamlit as st #pip install streamlit
import shutil
import random
import smtplib, ssl
from email.message import EmailMessage
import mimetypes
import paramiko
import os
from genconfigfiles import *
from resources import *
import logging

def sendEmail(message,session_dict,attach_label:bool):
    msg = EmailMessage()
    # smtp_server = "smtp.gmail.com"
    receiver_email = session_dict['user_email']+"@gmail.com"
    msg['From']=sender_email
    msg['To']=receiver_email
    msg['Subject']="The Sipros message"
    msg.set_content(message)
    if attach_label:
        filename = os.path.join(session_dict["whole_path"],"results.tar.gz")
        mime_type, _ = mimetypes.guess_type(filename)
        mime_type, mime_subtype = mime_type.split('/', 1)
        with open(filename, 'rb') as ap:
            msg.add_attachment(ap.read(), maintype=mime_type, subtype=mime_subtype,filename="results.tar.gz")
        ap.close()
    context = ssl.create_default_context()
    server = smtplib.SMTP(smtp_server,port)
    server.ehlo() # Can be omitted
    server.starttls(context=context) # Secure the connection
    server.ehlo() # Can be omitted
    server.login(sender_email, password)
    server.send_message(msg)
    server.quit() 

def Reversedatabase(inputfile:str, outputfile:str,config_dict:dict):
    probability_1 = 0.5
    probability_2 = 1
    
    training_prefix_str = 'Rev1_'
    testing_prefix_str = 'Rev_'
    
    if config_dict["Search_Type"] == 'SIP': # no training is needed for SIP search, so only test decoy is generated
        probability_1 = -1
        
    
    output_file = open(outputfile, "w")
    id_str = ""
    seq_str = ""
    seq_new_str = ""
    input_file = open(inputfile, "r")
    line_str = ""
    random_float = 0.0
    for line_str in input_file:
        if line_str[0] == '>':
            if seq_str != "":
                seq_new_str = (seq_str[::-1])
                output_file.write(id_str)
                output_file.write(seq_str)
                output_file.write('\n')
                random_float = random.random()
                if random_float <= probability_1:
                    output_file.write('>')
                    output_file.write(training_prefix_str)
                elif random_float <= probability_2:
                    output_file.write('>')
                    output_file.write(testing_prefix_str)                
                output_file.write(id_str[1:])
                output_file.write(seq_new_str)
                output_file.write("\n")
            id_str = line_str
            seq_str = ""
        else:
            seq_str += line_str.strip()

    if seq_str != "":
        seq_new_str = (seq_str[::-1])
        output_file.write(id_str)
        output_file.write(seq_str)
        output_file.write('\n')
        random_float = random.random()
        if random_float <= probability_1:
            output_file.write('>')
            output_file.write(training_prefix_str)
        elif random_float <= probability_2:
            output_file.write('>')
            output_file.write(testing_prefix_str)                 
        output_file.write(id_str[1:])
        output_file.write(seq_new_str)
        output_file.write("\n")
    id_str = line_str
    seq_str = ""
        
    input_file.close()
    output_file.close()

def Checkfiles(files_list:list):
    label_data = False
    count_fasta = 0
    for one_file in files_list:
        if one_file.find(".fa")!=-1 or one_file.find(".fasta")!=-1:
            count_fasta += 1
            print(one_file)
        # elif one_file.find(".ms2")!=-1 or one_file.find(".ft2")!=-1 or one_file.find(".mzML")!=-1 or one_file.find(".raw")!=-1: #############################
        if one_file.find(".raw") != -1 or one_file.find(".RAW") !=-1:
            label_data = True
    if label_data and count_fasta == 1:
        return True
    else:
        return False

def Generateconfigfile(session_dict:dict,times:int, parallel_label,search_type):
    data_file2_size = 0
    data_file1_size = 0
    fasta_size = 0
    ms_file = ""
    for tmp in os.listdir(session_dict['whole_path']):
        # if tmp.find(".ft2")!=-1 or tmp.find(".ms2") != -1 or tmp.find(".mzML")!=-1:
        if tmp.find(".FT2") !=-1 :
            data_file2_size += os.path.getsize(os.path.join(session_dict['whole_path'],tmp))
            ms_file = tmp
        elif tmp.find("FT1")!=-1:
            data_file1_size += os.path.getsize(os.path.join(session_dict['whole_path'],tmp))
        elif tmp.find(".fa")!=-1:
            fasta_size = os.path.getsize(os.path.join(session_dict['whole_path'],tmp))
    data_file2_size = (data_file2_size/1048576)/177
    data_file1_size = (data_file1_size/1048576)/106
    fasta_size = (fasta_size/1048576)/188
    if data_file2_size < 1:
        data_file2_size=1
    if fasta_size < 1:
        fasta_size = 1
    if data_file1_size<1:
        data_file1_size=1
    total_time = (data_file2_size*fasta_size*10*1.2+data_file1_size*1*1.2)*times

    path =os.path.join(session_dict['whole_path'],"sbatchscript.sh")
    if parallel_label:
        total_time = total_time/8
    hour = int(total_time/60)
    mins= int(total_time%60)
    if hour == 0 and mins<30:
        mins=30
    elif hour >= 48:
        hour = 47
        mins = 59
    params={"hour":hour, "mins":mins, "ms file":ms_file, "search type":search_type}

    if not parallel_label:
        generateSlurmOpenmp(path,params,session_dict)
    else:
        path2 =os.path.join(session_dict['whole_path'],"sbatchscript2.sh")
        generateSlrumMPI(path,path2,params,session_dict)
        
    return parallel_label

def generateUploadFile(path:str):
    file="label.txt"
    with open(os.path.join(path,file),"w") as f:
        f.write("loading:0\n")
        f.write("multiple:0\n")
        f.write("files:\n")
    f.close()

def writeUploadFile(path:str,content:tuple,files:str):
    file="label.txt"
    result=""
    with open(os.path.join(path,file),"r") as f:
        for line in f:
            if line.split(":")[0]==content[0]:
                result+=content[0]+":"+content[1]+"\n"
            elif len(files)!=0 and line.split(":")[0]=="files":
                result+="files:"+files+"\n"
            else:
                result+=line

    f.close()
    f=open(os.path.join(path,file),"w")
    f.write(result)
    f.close()

def readUploadFile(path:str,label:str):
    file="label.txt"
    result=""
    with open(os.path.join(path,file),"r") as f:
        for line in f:
            tmp=line.strip().split(":")
            if tmp[0]==label:
                result=tmp[1]
                break
    f.close()
    return result



def Generateuserinfo(max_jobs:int,path:str,user_email:str):
    file = "config.txt"
    with open(os.path.join(path,file),"w") as f:
        for i in range(max_jobs):
            f.write(str(i)+"\t0\t0\n")
        f.write(user_email)
    f.close()

def Readuserinfo(path:str):
    user_info_dict = {}
    file = "config.txt"
    with open(os.path.join(path,file),"r") as f:
        for line in f:
            tmp = line.strip().split("\t")
            user_info_dict[tmp[0]]=[tmp[1],tmp[2]]
            break
    f.close()
    return user_info_dict

def Writeuserinfo(session_dict,path:str):
    file = "config.txt"
    with open(os.path.join(path,file),"w") as f:
        for id in session_dict['user_info'].keys():
            f.write(str(id)+"\t"+str(session_dict['user_info'][id][0])+"\t"+str(session_dict['user_info'][id][1])+"\n")
        f.write(session_dict["user_email"])
    f.close()

def Checkavab(user_info)->str:
    for id in user_info.keys():
        if user_info[id][0] == "0":
            return id
    return "-1"

def Runonserver(host_inf:dict,session_dict:dict,del_label:bool,sbatch_config_file:str, again_label:bool):
    logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.INFO)
    multiple_label = readUploadFile(os.path.join(user_data_path,session_dict['user_email']),"multiple")
    client = paramiko.client.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(host_inf["host"], username=host_inf["username"], password=host_inf["passwords"])
    # _,current_path,_ = client.exec_command("pwd")
    # current_path = current_path.read().decode().strip()
    the_path = os.path.join("/","scratch","jiancheng",session_dict['whole_path'].split("/")[-2]+"+"+session_dict['whole_path'].split("/")[-1])
    if again_label:
        if not del_label:
            _,bb,bbb=client.exec_command("mkdir "+the_path)
            print(bbb.read().decode().strip())
            _,bb,bbb=client.exec_command("chmod -R 777 "+the_path)
            # print(bbb.read().decode().strip())
            if multiple_label=="1":
                _,tmp_path,tmp_err=client.exec_command("mkdir "+os.path.join(the_path,"raw_data"))
                # _,bb,bbb=client.exec_command("chmod -R 777 "+os.path.join(the_path,"raw_data"))
                # print("create folder for raw data")
                print(tmp_err.read().decode().strip())
                # print(tmp_path.read().decode().strip())
                # print(tmp_err.read().decode().strip())
        else:
            client.exec_command("rm -r"+os.path.join(the_path,"result"))
        # print("result")
        _,bb,bbb=client.exec_command("mkdir "+os.path.join(the_path,"result"))
        print(bbb.read().decode().strip())
        _,bb,bbb=client.exec_command("chmod -R 777 "+os.path.join(the_path,"result"))
        # print(bbb.read().decode().strip())
        sftp_client=client.open_sftp()
        if not del_label:
            for one_file in os.listdir(session_dict['whole_path']):
                os.system("chmod 777 "+os.path.join(session_dict['whole_path'],one_file))
                if multiple_label=="1" and one_file.find(".FT2") != -1:
                    _,bb,bbb=client.exec_command("chmod -R 777 "+os.path.join(the_path,"raw_data"))
                    # print(bbb.read().decode().strip())
                    sftp_client.put(os.path.join(session_dict['whole_path'],one_file),os.path.join(the_path,"raw_data",one_file))
                elif one_file.find(".FT1") != -1:
                    if multiple_label=="1":
                        _,aa,aaa=client.exec_command("mkdir "+os.path.join(the_path,one_file[:-4]))
                        print(aaa.read().decode().strip())
                        _,aa,aaa=client.exec_command("chmod -R 777 "+os.path.join(the_path,one_file[:-4]))
                        # print(aa.read().decode().strip())
                        # print(aaa.read().decode().strip())
                        try:
                            sftp_client.put(os.path.join(session_dict['whole_path'],one_file),os.path.join(the_path,one_file[:-4],one_file))
                        except:
                            print("error")
                            print(os.path.join(session_dict['whole_path'],one_file))
                            print(os.path.join(the_path,one_file[:-4],one_file))
                    else:
                        sftp_client.put(os.path.join(session_dict['whole_path'],one_file),os.path.join(the_path,"result",one_file))
                else:
                    sftp_client.put(os.path.join(session_dict['whole_path'],one_file),os.path.join(the_path,one_file))
        sftp_client.close()
    _, aa, aaa = client.exec_command("ls "+the_path)
    print("super:")
    print(aa.read().decode().strip())
    in_job, job_id, error_job = client.exec_command("cd "+the_path+"; sbatch "+sbatch_config_file)
    job_id = job_id.read().decode()
    error_job = error_job.read().decode()
    logging.error(error_job)
    logging.info("job information:"+session_dict['user_email'] + ":" +job_id)
    session_dict['job_id'] = job_id.split(" ")[3].strip()
    
    client.close()  


def RunDemoServer(host_inf:dict,local_path:str,del_label:bool,sbatch_config_file:str, again_label:bool)->str:
    logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.INFO)
    # multiple_label = readUploadFile(os.path.join(user_data_path,session_dict['user_email']),"multiple")
    client = paramiko.client.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(host_inf["host"], username=host_inf["username"], password=host_inf["passwords"])
    # _,current_path,_ = client.exec_command("pwd")
    # current_path = current_path.read().decode().strip()
    the_path = os.path.join("/","scratch","jiancheng","example")
    if again_label:
        # if not del_label:
        _,bb,bbb=client.exec_command("mkdir "+the_path)
        print(bbb.read().decode().strip())
        _,bb,bbb=client.exec_command("chmod -R 777 "+the_path)
        print(bbb.read().decode().strip())
            # if multiple_label=="1":
            #     _,tmp_path,tmp_err=client.exec_command("mkdir "+os.path.join(the_path,"raw_data"))
            #     # _,bb,bbb=client.exec_command("chmod -R 777 "+os.path.join(the_path,"raw_data"))
            #     # print("create folder for raw data")
            #     print(tmp_err.read().decode().strip())
            #     # print(tmp_path.read().decode().strip())
            #     # print(tmp_err.read().decode().strip())
        # else:
            # client.exec_command("rm -r"+os.path.join(the_path,"result"))
        # print("result")
        _,bb,bbb=client.exec_command("mkdir "+os.path.join(the_path,"result"))
        print(bbb.read().decode().strip())
        _,bb,bbb=client.exec_command("chmod -R 777 "+os.path.join(the_path,"result"))
        # print(bbb.read().decode().strip())
        sftp_client=client.open_sftp()
        if not del_label:
            for one_file in os.listdir(local_path):
                os.system("chmod 777 "+os.path.join(local_path,one_file))
                # if multiple_label=="1" and one_file.find(".FT2") != -1:
                #     _,bb,bbb=client.exec_command("chmod -R 777 "+os.path.join(the_path,"raw_data"))
                #     # print(bbb.read().decode().strip())
                #     sftp_client.put(os.path.join(session_dict['whole_path'],one_file),os.path.join(the_path,"raw_data",one_file))
                if one_file.find(".FT1") != -1:
                    # if multiple_label=="1":
                    #     _,aa,aaa=client.exec_command("mkdir "+os.path.join(the_path,one_file[:-4]))
                    #     print(aaa.read().decode().strip())
                    #     _,aa,aaa=client.exec_command("chmod -R 777 "+os.path.join(the_path,one_file[:-4]))
                    #     # print(aa.read().decode().strip())
                    #     # print(aaa.read().decode().strip())
                    #     try:
                    #         sftp_client.put(os.path.join(session_dict['whole_path'],one_file),os.path.join(the_path,one_file[:-4],one_file))
                    #     except:
                    #         print("error")
                    #         print(os.path.join(session_dict['whole_path'],one_file))
                    #         print(os.path.join(the_path,one_file[:-4],one_file))
                    # else:
                    sftp_client.put(os.path.join(local_path,one_file),os.path.join(the_path,"result",one_file))
                else:
                    sftp_client.put(os.path.join(local_path,one_file),os.path.join(the_path,one_file))
        sftp_client.close()
    _, aa, aaa = client.exec_command("ls "+the_path)
    print("super:")
    print(aa.read().decode().strip())
    in_job, job_id, error_job = client.exec_command("cd "+the_path+"; sbatch "+sbatch_config_file)
    job_id = job_id.read().decode()
    error_job = error_job.read().decode()
    logging.error(error_job)
    logging.info("job information: example -" +job_id)
    
    client.close() 
    return job_id.split(" ")[3].strip()

def copyFiles(host_inf:dict,session_dict:dict,search_type:str)->None:
    file=readUploadFile(os.path.join(user_data_path,session_dict['user_email']),"files").strip().split("~")
    client = paramiko.client.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(host_inf["host"], username=host_inf["username"], password=host_inf["passwords"])
    # _,current_path,_ = client.exec_command("pwd")
    # current_path = current_path.read().decode().strip()
    the_path = os.path.join("/","scratch","jiancheng",session_dict['whole_path'].split("/")[-2]+"+"+session_dict['whole_path'].split("/")[-1])
    for tmp in file:
        if search_type.find("SIP")!=-1:
            client.exec_command("cp "+the_path+"/result/*.pro2psm.cluster.txt "+the_path+"/"+tmp+"/"+tmp+".pro2psm.txt")
        else:
            client.exec_command("cp "+the_path+"/result/*.pro2psm.txt "+the_path+"/"+tmp+"/"+tmp+".pro2psm.txt")
    client.close()


def Checkstatus(host_inf:dict, job_id:int,whole_path:str):
    client = paramiko.client.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(host_inf["host"], username=host_inf["username"], password=host_inf["passwords"])
    _,status,_ = client.exec_command("sacct --format=jobid,elapsed,state --jobs="+job_id)
    status = status.read().decode().split("\n")

    _,files,_ = client.exec_command("ls "+whole_path.split("/")[-2]+"+"+whole_path.split("/")[-1])
    client.close()
    for one_line in status[2:]:
        status_list = []
        tmp_word = ""
        start_label = True
        for c in one_line:
            if c != " ":
                tmp_word += c
                start_label = True
            elif start_label and c ==" ":
                status_list.append(tmp_word)
                tmp_word = ""
                start_label = False
        if status_list[0] == job_id:
            if status_list[2]=="COMPLETED":
                for one_file in files:
                    if one_file.find("_stderr.txt") != -1 or one_file.find("_stdout.txt")!=-1:
                        with open(one_file,"r") as f:
                            for line in f:
                                if line.find("Program exit!")!=-1 or line.find("ERROR!")!=-1:
                                    f.close()
                                    return "FAILED", status_list[1]
                        f.close()
            return status_list[2], status_list[1]

        return "",""

def Downloadfile(host_inf:dict,session_dict:dict,search_type:str):
    files=readUploadFile(os.path.join(user_data_path,session_dict['user_email']),"files").strip().split("~")

    client = paramiko.client.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(host_inf["host"], username=host_inf["username"], password=host_inf["passwords"])
    sftp_client = client.open_sftp()
    # _,current_path,_ = client.exec_command("pwd")
    # current_path = current_path.read().decode().strip()
    the_path = os.path.join("/","scratch","jiancheng",session_dict['whole_path'].split("/")[-2]+"+"+session_dict['whole_path'].split("/")[-1])
    for tmp in files:
        _,tar_state,_=client.exec_command("cp "+the_path+"/"+tmp+"/"+"*Peptide.txt "+the_path+"/result")
        _,tar_state,_=client.exec_command("cp "+the_path+"/"+tmp+"/"+"*Protein.txt "+the_path+"/result")
    if search_type == "SIP":
        _,tar_state,_=client.exec_command("cd "+the_path+"/result"+" ; "+"tar -cvzf results.tar.gz *.txt ; chmod 777 results.tar.gz")

    else:
        _,tar_state,bbb=client.exec_command("cd "+the_path+"/result"+" ; "+"tar -cvzf results.tar.gz *.txt ; chmod 777 results.tar.gz")
        print(bbb.read().decode().strip())

    # print(tar_state.read().decode().strip())
    sftp_client.get(os.path.join(the_path,"result","results.tar.gz"),os.path.join(session_dict['whole_path'],"results.tar.gz"))
    sftp_client.close()
    client.close()

def downloadDemoFile(host_inf:dict):
    client = paramiko.client.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(host_inf["host"], username=host_inf["username"], password=host_inf["passwords"])
    sftp_client = client.open_sftp()
    # _,current_path,_ = client.exec_command("pwd")
    # current_path = current_path.read().decode().strip()
    the_path = os.path.join("/","scratch","jiancheng","example")
    _,files,_ = client.exec_command("ls "+the_path+"/result")
    for tmp in files:
        sftp_client.get(os.path.join(the_path,"result",tmp),os.path.join("/home/jiancheng/example",tmp))
    sftp_client.close()
    client.close()

def delOnSlurm(host_inf:dict, session_dict:dict):
    client = paramiko.client.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(host_inf["host"], username=host_inf["username"], password=host_inf["passwords"])
    the_path = os.path.join("/","scratch","jiancheng",session_dict['whole_path'].split("/")[-2]+"+"+session_dict['whole_path'].split("/")[-1])
    _,status,_ = client.exec_command("rm -r "+the_path)
    client.close()
    shutil.rmtree(session_dict['whole_path'])
