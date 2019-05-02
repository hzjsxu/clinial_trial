# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 13:51:57 2019

@author: jsxu
"""

from __future__ import division
from collections import defaultdict
import xml.etree.ElementTree as ET
import re, json, nltk

def process_node(parent, parent_dict):
    if list(parent):
        for child in parent:
            child_dict = {}
            child_dict = process_node(child, child_dict)
            parent_dict[child.tag] = child_dict
    else:
        parent_dict['field_value'] = parent.text
    return parent_dict

def load_trial(trial_id, root_dir = 'E:/bio_TM/clinical_trials_data/'):
    dir = trial_id[:-4] + 'xxxx'
    filepath = root_dir + dir + '/' + trial_id + '.xml'
    tree = ET.parse(filepath)
    root = tree.getroot()
    root_dict = {}
    return process_node(root, root_dict)

def find_textblocks(field, textblocks = [], current = 'root'):
    subfields = field.keys()
    if not subfields == ['field_value']:
        if 'textblock' in subfields:
            textblocks.append((current, field['textblock']['field_value']))
        if len(subfields) > 1:
            for subfield in subfields:
                if not subfield == 'textblock':
                    textblocks = find_textblocks(field[subfield], textblocks, current = current + '>' + subfield)
    return textblocks

def process_textblock(textblock):
    utext = textblock
    lines = re.sub('\n*[ ]+', ' ', utext)
    lines = lines.strip()
    return lines
    
def get_textblocks(textblocks):
    textblocks_dict = {}
    for textblock in textblocks:
        levels = re.split('>', textblock[0])
        key = '>'.join(levels[1:])
        textblocks_dict[key] = process_textblock(textblock[1])
    return textblocks_dict

def extract_textblocks(trial_id, root_dir = 'E:/bio_TM/clinical_trials_data/'):
    return get_textblocks(find_textblocks(load_trial(trial_id, root_dir = root_dir), textblocks = []))

def create_regex_format(x):
    regex = re.compile(r"[^a-z0-9](" + x + r")[^a-z0-9]")
    return regex
    
def find_indicators(element, indicators_dict = {}, current = "root"):
    if type(element) == dict:
        subelements = element.keys()
        if "indicators" in subelements and element["indicators"]:
            indicators_dict[re.sub("refinements>", "",current)] = (map(create_regex_format, element["indicators"]))
        if len(subelements) > 1:
            for subelement in subelements:
                if type(element[subelement]) == dict:
                    indicators_dict = find_indicators(
                        element[subelement], indicators_dict, current = current + ">" + subelement
                    )
    return indicators_dict

def match_indicators(d, indicators_dict, secID):
    lines = d.get(secID,"")
    found = defaultdict(list)
    for sentNum, sentence in enumerate(nltk.sent_tokenize(lines)):
        for group in indicators_dict:
            matches = []
            for n, indicator in enumerate(indicators_dict[group]):
                padded = " " + sentence.lower().strip() + " "
                if indicator.search(padded):
                    nuggets = indicator.findall(padded)
                    matches.append(re.sub(r"^.*?\]\((.*?)\)\[.*?$",r"\1",indicator.pattern))
            if matches:
                found[group].append((matches, str(sentNum)))
    return dict(found)


schema = json.load(open("E:/bio_TM/structured-output-schema.json", "r"))
indicators_dict = find_indicators(schema)
textblocks_dict = extract_textblocks("NCT02775682")
found = match_indicators(textblocks_dict, indicators_dict, "brief_summary")
print(found)
















