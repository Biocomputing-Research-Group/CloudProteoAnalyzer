# CloudProteoAnalyzer
CloudProteoAnalyzer is developed to provide a convenient online service for proteomics data analysis, including identification and quantification. It is implemented within a high-performance computing (HPC) cluster to enhance computational scalability for database searching, which is distributed across multiple computing nodes and multi-threaded on multi-core CPUs. The user-friendly web interface of CloudProteoAnalyzer empowers users to customize all parameters for protein identification and quantification.

The CloudProteoAnalyzer link is https://sipros.oscer.ou.edu/ or https://sipros-oscer.unt.edu.
## Folder Description 
### CloudProteoAnalyzer Folder
There is the source code of the CloudProteoAnalyzer
### Script Folder
There are script files for obtaining benchmark results in the paper.
### Configuration Folder
The parameters of the methods are used for benchmark results in the paper.
### Database Folder
The protein databases are used for benchmark results in the paper.
## Software Setup
### Dependent Software
* Apache 2.4.41
* Raxport (https://github.com/xyz1396/Raxport.net)
* Slurm (https://slurm.schedmd.com/documentation.html)
* Sipros Ensemble (https://github.com/thepanlab/Sipros4)
* ProRata (Download from the PraRata folder)
### Dependencies
* Python 3.8.10
* Python libraries: Streamlit, Numpy, Pandas, Scipy, scikit-learn, google-api-python-client, google-auth-httplib2, google-auth-oauthlib, Requests, SSL, and Paramiko.
## Installation
### Supercomputer Part
1. Download and install Slurm
2. Download and install Sipros Ensemble
3. Download ProRata
4. Add the path of the executing files of Sipros Ensemble and ProRata to the system path
### Web Application Server Part
1. Download and install Apache
2. Download and install Raxport
3. Download the source code
4. Create the Google client ID credentials (https://developers.google.com/workspace/guides/create-credentials) and download the credential JSON file
5. Setup some configurations in resoures.py
6. Run server
   ```
   streamlit run ui_db.py --server.address localhost --server.port 8502 --server.enableXsrfProtection=false --server.enableCORS=false . &
   ```
