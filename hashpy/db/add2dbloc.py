#!/usr/bin/env python
#%run dbhash.py -r --pf dbhash.pf --evid 383991  /home/markwill/working/dbhash/reno
#
from hashpy.db.utils import add_antelope_path
add_antelope_path()
from antelope.stock import *
import sys,os

f = os.path.join(os.environ['ANTELOPE'],'data','pf','dbloc2.pf')

# note: pfget('dbloc2', 'User') will also work, but will traverse pfpath
new_User = pfget(f, 'User')
new_Env  = pfget(f, 'Env')

omi = new_User['origin_menu_items']

t = tuple([x for x in omi if 'dbhash' not in x] + ['Calculate_focal_mech dbhash --loc'])
new_User['origin_menu_items'] = t

if 'PRESERVE' not in new_Env['PATH']:
	new_Env['PATH'] = 'PRESERVE ||' + new_Env['PATH']

pfput('User', new_User, f)
pfput('Env', new_Env, f)
