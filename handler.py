#!/usr/bin/env python
# Copyright 2016 Station X, Inc. or its affiliates. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import ctypes
import json
import os

# use python logging module to log to CloudWatch
# http://docs.aws.amazon.com/lambda/latest/dg/python-logging.html
import logging
logging.getLogger().setLevel(logging.DEBUG)

# must load all shared libraries and set the R environment variables before we can import rpy2
# load R shared libraries from lib dir
for file in os.listdir('lib'):
    if os.path.isfile(os.path.join('lib', file)):
        ctypes.cdll.LoadLibrary(os.path.join('lib', file))

# set R environment variables
os.environ["R_HOME"] = os.getcwd()
os.environ["R_LIBS"] = os.path.join(os.getcwd(), 'site-library')

# import rpy2
import rpy2
from rpy2 import robjects
from rpy2.robjects import r

def calculate_survival_stats(times, events, values_by_record):
    
    logging.debug('Calculating stats')

    #load R library
    r.source("script.R")

    logging.debug('Done calculating stats')

def lambda_handler(event, context):
    times = event['times']
    events = event['events']
    # support receiving values (ex: expression) for multiple records (ex: genes)
    values_by_record = event['values_by_record']
    logging.info('Number of samples: {0}'.format(len(times)))
    logging.info('Number of genes/variants: {0}'.format(len(values_by_record)))

    try:
        stats_list = calculate_survival_stats(times, events, values_by_record)
    except rpy2.rinterface.RRuntimeError, e:
        logging.error('Payload: {0}'.format(event))
        logging.error('Error: {0}'.format(e.message))

        # generate a JSON error response that API Gateway will parse and associate with a HTTP Status Code
        error = {}
        error['errorType'] = 'StatisticsError'
        error['httpStatus'] = 400
        error['request_id'] = context.aws_request_id
        error['message'] = e.message.replace('\n', ' ') # convert multi-line message into single line
        raise Exception(json.dumps(error))

    res = {}
    res['statistics_list'] = stats_list
    return res
