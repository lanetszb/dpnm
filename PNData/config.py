import sys
import os
import configparser

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))

config = configparser.ConfigParser()
config.read(filenames=sys.argv[1])

for section in config.sections():
    print()
    print('[' + section + ']')
    for key in config[section]:
        print(key, config.get(section, key))
