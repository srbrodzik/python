#!/usr/bin/python3

import requests

def valid(url):
    get = requests.get(url)
    if get.status_code == 200:
        return True
    else:
        return False
