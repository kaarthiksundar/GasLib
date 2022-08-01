import argparse
from logger import CustomFormatter
import logging
from fileprocessor import parse_gaslib


def setup_argparse(): 
    parser = argparse.ArgumentParser()
    parser.add_argument('--datafolder', help='folder where GasLib data is available', default='./data/')
    parser.add_argument('--file', help='GasLib zipped directory', default='GasLib-11.zip')
    return parser

if __name__ == '__main__':
    # create logger with 'spam_application'
    log = logging.getLogger("GasLib")
    log.setLevel(logging.DEBUG)

    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(CustomFormatter())
    
    log.addHandler(ch)

    args = setup_argparse().parse_args()
    datafile = args.datafolder + args.file 
