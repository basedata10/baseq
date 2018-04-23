import os, sys
import basematic.setting as setting
import configparser as CP

class configManager:
    def __init__(self, path = setting.Config):
        if not os.path.exists(path):
            sys.exit("[error] Config File Do Not Exist, please run: basematic init")
        self.file = path
        self.config = CP.ConfigParser()
        self.config.optionxform = str
        self.config.read(path)

    def get_section(self, section):
        res = {}
        if section not in self.config:
            sys.exit("[error] Section [{}] not defined in config file '{}'".format(section, self.file))
        for key in self.config[section]:
            res[key] = self.config[section][key]
        return sectionManager(res, section, self)

class sectionManager:
    def __init__(self, data, name, config):
        self.data = data
        self.name = name
        self.keys = self.get_keys()
        self.config = config

    def get_keys(self):
        keys = []
        for key in self.data:
            keys.append(key)
        return keys

    def get(self, key):
        if key in self.keys:
            print("[info] [{}] {} : {}".format(self.name, key, self.data[key]))
            return self.data[key]
        else:
            sys.exit("[error] '{}' is not configured in [{}]. \n[error] Please add the item in '{}'.".format(key, self.name, self.config.file))

def get_config(section, item):
    section_cfg = configManager().get_section(section)
    cfg = section_cfg.get(item)
    return cfg