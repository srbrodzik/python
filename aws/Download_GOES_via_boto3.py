#!/usr/bin/python3

import requests
import boto3
from botocore import UNSIGNED
from botocore.client import Config

def get_s3_keys(bucket, s3_client, prefix = ''):
    """
    Generate the keys in an S3 bucket.
    :param bucket: Name of the S3 bucket.
    :param prefix: Only fetch keys that start with this prefix (optional).
    """
    kwargs = {'Bucket': bucket}

    if isinstance(prefix, str):
        kwargs['Prefix'] = prefix

    while True:
        resp = s3_client.list_objects_v2(**kwargs)
        for obj in resp['Contents']:
            key = obj['Key']
            if key.startswith(prefix):
                yield key

        try:
            kwargs['ContinuationToken'] = resp['NextContinuationToken']
        except KeyError:
            break

# Inputs
bucket_name = 'noaa-goes16'
product_name = 'ABI-L1b-RadF'
year = 2019
day_of_year = 79
hour = 14
band = 3
download_file = True

# Initialize s3 client
s3_client = boto3.client('s3', config=Config(signature_version=UNSIGNED))

# Get keys (full path names) of files meeting input criteria
keys = get_s3_keys(bucket_name,
                   s3_client,
                   prefix = f'{product_name}/{year}/{day_of_year:03.0f}/{hour:02.0f}/OR_{product_name}-M3C{band:02.0f}')

# Select the first measurement taken within the hour​ (as an example)​
key = [key for key in keys][0]
file_name = key.split('/')[-1].split('.')[0]

resp = requests.get(f'https://{bucket_name}.s3.amazonaws.com/{key}')
if download_file:
    open(file_name,'wb').write(resp.content)
