import boto3
from datetime import datetime, timedelta, tzinfo

def lambda_handler(event, context):
    # AWS Batch job and S3 bucket definitions:
    JOB_QUEUE="spot40_default"
    JOB_DEFINITION="birdcast-observed-images:1"
    BUCKET="vol2bird"
    S3_PROFILES="output"
    S3_IMAGES="mosaic"
    DATETIME = datetime.utcnow() - timedelta(minutes = 15)

    # construct the AWS Batch parameters to parse
    job_params={'date':'--date='+DATETIME.strftime('%Y%m%d%H%M'), 'bucket':'--bucket='+BUCKET,'s3_profiles':'--s3_profiles='+S3_PROFILES,'s3_images':'--s3_images='+S3_IMAGES, 'log_filenames':'--log_filenames'}
    # name the job
    job_name=DATETIME.strftime('mosaic_%Y%m%d%H%M')
    
    # connect to Batch
    client = boto3.client('batch')
    
    # submit the job
    response = client.submit_job(
       jobDefinition=JOB_DEFINITION,
       jobName=job_name,
       jobQueue=JOB_QUEUE,
       parameters=job_params)
