##
# alignment_server.py
# @created 18 April 2016
# @author Kyle R. Covington; @email kylecovington1@gmail.com
#
# provides an http server to access multiz alignment conversions
# once running, users can do something like 
# @code curl -d '{"source": {"species": "hg19", "chrom": "1", "start": 100, "end": 110}, "target": "mm10"}' http://localhost:6480
# 

import os, os.path, sys
import json, glob, threading, time, gzip
import alignment
from BaseHTTPServer import BaseHTTPRequestHandler, HTTPServer
import SocketServer
import logging

logger = logging.getLogger('malign.alignment_server')
logger.addHandler(logging.NullHandler())
logger.setLevel(logging.DEBUG)


class AlignmentHelper(object):
    def __init__(self, maf_dir, subject = 'hg19', query = 'mm10'):
        self.alignment_dict = {}
        self.subject = subject
        self.query = query
        self.alignment_dict_lock = threading.RLock()
        self.available_alignments = {os.path.basename(x).replace('chr', '').split('.')[0]: x for x in glob.glob(os.path.join(maf_dir, '*.maf.gz'))}
    def load_alignment(self, chrom):
        if chrom in self.alignment_dict:
            return
        elif chrom in self.available_alignments:
            with self.alignment_dict_lock:
                logger.info("Loading chromosome %s", chrom)
                with gzip.open(self.available_alignments[chrom], 'rb') as fi:
                    self.alignment_dict[chrom] = alignment.AlignmentContainer(fi, self.subject, self.query)
                logger.info("Chrom %s loaded", chrom)
        else:
            raise ValueError("chrom is not valid, perform availablechroms to get available chromosomes")
    def availablechroms(self):
        return self.available_alignments.keys()

    def get_alignment(self, chrom, start, end):
        if chrom in self.alignment_dict:
            return self.alignment_dict[chrom].get_alignment(start, end)
        else:
            raise ValueError("chrom %s has not been loaded" % chrom)



class S(BaseHTTPRequestHandler):
    helper = None
    def _set_headers(self):
        self.send_response(200)
        self.send_header("Content-type", "application/json")
        self.end_headers()

    def do_HEAD(self):
        self._set_headers()

    def do_POST(self):
        self._set_headers()
        length = int(self.headers.getheader('content-length'))
        input = json.loads(self.rfile.read(length))
        # #####
        # now we have to decide how to process the request
        # #####
        try:
            if self.path == '/availablechroms':
                response = self._get_available_chroms()
            elif self.path == '/sequence':
                response = self._get_sequence(input)
            else:
                raise ValueError("path %s is not known" % self.path)
        except Exception as inst:
            response = {'action': 'rejected', 
                        'reason': str(inst), 
                        'options': [
                            '/availablechroms - lists available chromosomes', 
                            '/sequence {"chrom": <>, "start": <>, "end": <>} - Returns the sequences at chrom start and end']}
        self.wfile.write(json.dumps(response))

    def _get_available_chroms(self):
        return {'action': 'accepted', 'chromosomes': self.helper.availablechroms()}

    def _get_sequence(self, input):
        for k in ('chrom', 'start', 'end'):
            if k not in input:
                raise ValueError("input is missing required key %s" % k)
        self.helper.load_alignment(input['chrom'])
        a = self.helper.get_alignment(input['chrom'], input['start'], input['end'])
        return {'action': 'accepted', 'sourceseq': a.get_source_sequence(input['start'], input['end']), 'queryseq': a.get_query_sequence(input['start'], input['end'])}

def main(args, server_class=HTTPServer, handler_class=S):
    server_address = ('', args.port)
    handler_class.helper = AlignmentHelper(args.MAFDIR, subject = args.subject, query = args.query)
    httpd = server_class(server_address, handler_class)
    logger.info("Server Started...")
    httpd.serve_forever()

if __name__ == '__main__':
    import argparse
    
    ch = logging.StreamHandler(stream=sys.stdout)
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    parser = argparse.ArgumentParser()
    
    parser.add_argument('--subject', type = str, default = 'hg19', help = 'subject id')
    parser.add_argument('--query', type = str, default = 'mm10', help = 'query id')
    parser.add_argument('--port', type = int, help = 'server port, default 6480', default = 6480)
    parser.add_argument('MAFDIR', type = str, help = 'directory with maf.gz files for processing')

    args = parser.parse_args()

    main(args)
