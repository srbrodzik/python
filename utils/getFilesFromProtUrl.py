from bs4 import BeautifulSoup
import requests
from requests.auth import HTTPBasicAuth

def get_flist(url,user,pwd,suffix):

    auth = HTTPBasicAuth(user,pwd)
    r = requests.get(url = url, auth=auth, verify=False)
    soup = BeautifulSoup(r.text, 'html.parser')
    flist = [url + '/' + node.get('href') for node in soup.find_all('a') if node.get('href').endswith(suffix)]
    
    return flist