'''
A simple authenticated web server handler
From: https://github.com/tianhuil/SimpleHTTPAuthServer
'''

from __future__ import print_function

import argparse
import base64
import os
import ssl
import sys

try:
    from SimpleHTTPServer import SimpleHTTPRequestHandler
except ImportError:
    from http.server import SimpleHTTPRequestHandler

try:
    from SocketServer import TCPServer
except ImportError:
    from socketserver import TCPServer

__prog__ = 'SimpleHTTPAuthServer'
__version__ = '1.3'

CERT_FILE = os.path.expanduser("~/.ssh/cert.pem")
KEY_FILE = os.path.expanduser("~/.ssh/key.pem")
SSL_CMD = "openssl req -newkey rsa:2048 -new -nodes -x509 "\
            "-days 3650 -keyout {0} -out {1}".format(KEY_FILE, CERT_FILE)

class SimpleHTTPAuthHandler(SimpleHTTPRequestHandler):
    ''' Main class to present webpages and authentication. '''
    KEY = ''

    def do_HEAD(self):
        ''' head method '''
        print("send header")
        self.send_response(200)
        self.send_header('Content-type', 'text/html')
        self.end_headers()

    def do_authhead(self):
        ''' do authentication '''
        print("send header")
        self.send_response(401)
        self.send_header('WWW-Authenticate', 'Basic realm=\"Test\"')
        self.send_header('Content-type', 'text/html')
        self.end_headers()

    def do_GET(self):
        ''' Present frontpage with user authentication. '''
        if self.headers.getheader('Authorization') is None:
            self.do_authhead()
            self.wfile.write('no auth header received')
        elif self.headers.getheader('Authorization') == 'Basic '+ self.KEY:
            SimpleHTTPRequestHandler.do_GET(self)
        else:
            self.do_authhead()
            self.wfile.write(self.headers.getheader('Authorization'))
            self.wfile.write('not authenticated')

def serve_https(https_port=80, https=True, start_dir=None, handler_class=SimpleHTTPAuthHandler):
    ''' setting up server '''
    httpd = TCPServer(("", https_port), handler_class)

    if https:
        httpd.socket = ssl.wrap_socket(httpd.socket, keyfile=KEY_FILE,
                                       certfile=CERT_FILE, server_side=True)

    if start_dir:
        print("Changing dir to {cd}".format(cd=start_dir))
        os.chdir(start_dir)

    socket_addr = httpd.socket.getsockname()
    print("Serving HTTP on", socket_addr[0], "port", socket_addr[1], "...")
    httpd.serve_forever()

def main():
    ''' Parsing inputs '''
    parser = argparse.ArgumentParser(prog=__prog__)
    parser.add_argument('port', type=int, help='port number')
    parser.add_argument('key', help='username:password')
    parser.add_argument('--dir', required=False, help='directory')
    parser.add_argument('--https', help='Use https', action='store_true', default=False)
    args = parser.parse_args()

    if args.https:
        if not (os.path.exists(CERT_FILE) and os.path.exists(KEY_FILE)):
            print("", file=sys.stderr)
            print("Missing {} or {}".format(CERT_FILE, KEY_FILE), file=sys.stderr)
            print("Run `{}`".format(SSL_CMD), file=sys.stderr)
            print("", file=sys.stderr)
            sys.exit(1)

    # print(args.key.encode('utf-8'))

    SimpleHTTPAuthHandler.KEY = base64.b64encode(args.key.encode('utf-8'))

    serve_https(int(args.port), https=args.https,
                start_dir=args.dir, handler_class=SimpleHTTPAuthHandler)


if __name__ == '__main__':
    main()