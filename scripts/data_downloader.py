#!/usr/bin/python

'''
Script name: data_downloader.py
Usage: python data_downloader.py (python2)
'''

import os


DATA_DIR = os.path.join("..", "data")
OUTPUT_FILE = os.path.join(DATA_DIR, "sequences")
URL = '"https://enhancer.lbl.gov/cgi-bin/imagedb3.pl?page_size=100;show=1;search.result=yes;search.form=no;form=search;action=search;search.org=Both;search.sequence=1"'

def download_data():
    if not os.path.exists(DATA_DIR):
        os.makedirs(DATA_DIR)
    if len(os.listdir(DATA_DIR)) == 0:
        os.system('wget -O {} {}'.format(OUTPUT_FILE, URL))

if __name__=="__main__":
    print "Running the script..."
    download_data()
    print "Done!"
