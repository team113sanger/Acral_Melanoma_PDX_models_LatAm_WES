# Ingestion and Initial QC using DERMATLAS ingestion pipeline



## Create a dermatlas_ingestion base:

```bash
PROJECTDIR=/lustre/scratch125/casm/teams/team113/projects/6633_2729_3248_PDX_models_from_Latin_America_WES
BASE_DIR=${PROJECTDIR}/base_dir

module load dermatlas-setup-ingestion-base

dermatlas-setup-ingestion-base ${PROJECTDIR}

# Staged location 2)
# /lustre/scratch125/casm/staging/

# Purge all of the modules from dermatlas-setup-ingestion-base
module purge

```

<!-- OPTIONAL NEXT STEPS

If you'd like to source-control your ingestion base, you can attach a remote
GitLab repository and push the changes:

cd /lustre/scratch125/casm/teams/team113/projects/6633_2729_3248_PDX_models_from_Latin_America_WES/base_dir/ingestion_base
git remote add origin <remote-repo-url>
git push -u origin develop
git add . '!scripts'
git commit -m 'Initial commit (from template)'
git push

Also consider ignoring the untracked files in the scripts submodule
and updating the README.md file. -->


## Next steps stage BAM files into Stage directory

```bash
PROJECTDIR=/lustre/scratch125/casm/teams/team113/projects/6633_2729_3248_PDX_models_from_Latin_America_WES
BASE_DIR=${PROJECTDIR}/base_dir/ingestion_base

cd ${PROJECTDIR}
source ${BASE_DIR}/source_me.sh

module load dermatlas-ingestion/0.0.7 
module load dataImportExport/1.58.2

STUDY=6633
PROJECT=2729
#STAGE_DIR="/lustre/scratch125/casm/staging/team113"
#STAGE_DIR="/lustre/scratch125/casm/staging/team113/pdx"
```

### Create a dir to stage into + a place to store the logs   

```bash	

mkdir -p "${STAGE_DIR}/${PROJECT}/logs"
```
### Get the sample list from the project
```bash
bash -c "mlwh-study-samples '$STUDY'; echo; " 2>/dev/null > "${STAGE_DIR}/${PROJECT}/sample_list.txt"
```
### Stage the BAM files into the stage directory
```bash
stageBam.pl -p ${PROJECT} \
-s "${STAGE_DIR}/${PROJECT}/sample_list.txt" \
-o "${STAGE_DIR}" \
-fo \
-lo "${STAGE_DIR}/${PROJECT}/logs" \
--type m
```

## Ingest the samples after BAM staging


Using `dermatlas-ingestion/0.0.7`

```bash
PROJECTDIR=/lustre/scratch125/casm/teams/team113/projects/6633_2729_3248_PDX_models_from_Latin_America_WES
BASE_DIR=${PROJECTDIR}/base_dir/ingestion_base
SRC_ME=${BASE_DIR}/source_me.sh
START_INGSCRIPT=${BASE_DIR}/scripts/start_ingestion.py
STUDY=6633
PROJECT=2729

source ${SRC_ME}

module load dermatlas-ingestion/0.0.7

#Check gogole drive authentication
t113-google-drive check-auth

#Check iRODS authentication and connection
seqscope irods iexpire

iinit

${START_INGSCRIPT} \
--source-me ${SRC_ME} \
--project ${PROJECT} \
--study ${STUDY} \
--dna

```
## NOTE: RNA ingestion requires to have a Google drive URL provided