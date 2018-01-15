#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Written by: Jimmy NGAI jimmycfngai@gmail.com

import logging

logging.basicConfig(level=logging.DEBUG,
        #format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%m-%d %H:%M:%S',
        filename='debug.log',
        filemode='a')   #mode append:a, overwritten:w 
# define a Handler which writes INFO messages or higher to the sys.stderr
console = logging.StreamHandler()
console.setLevel(logging.INFO)
# set a format which is simpler for console use
#formatter = logging.Formatter('%(levelname)-8s : %(message)s')
formatter = logging.Formatter('%(message)s')
# tell the handler to use this format
console.setFormatter(formatter)
logger = logging.getLogger(__name__)
# add the handler to the root logger
logger.addHandler(console)
#TODO: to change the logger level from cmd args.
#TODO: and disable debug.log file log
