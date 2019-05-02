# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 21:32:43 2019

@author: jsxu
"""

from __future__ import division
from collections import defaultdict
import xml.etree.ElementTree as ET
import re, json, nltk, random, spacy

nlp = spacy.load('en_core_web_sm')
def process_node(parent, parent_dict):
    if list(parent):
        for child in parent:
            child_dict = {}
            child_dict = process_node(child, child_dict)
            parent_dict[child.tag] = child_dict
    else:
        parent_dict["field_value"] = parent.text
    return parent_dict

def load_trial(trial_id, root_dir = "E:/bio_TM/clinical_trials_data/"):
    directory = trial_id[:-4] + "xxxx"
    filepath = root_dir + directory + "/" + trial_id + ".xml"
    tree = ET.parse(filepath)
    root = tree.getroot()
    root_dict = {}
    return process_node(root, root_dict)

def find_textblocks(field, textblocks = [], current = "root"):
    subfields = field.keys()
    if not subfields == ["field_value"]:
        if "textblock" in subfields:
            textblocks.append((current, field["textblock"]["field_value"]))
        if len(subfields) > 1:    
            for subfield in subfields:
                if not subfield == "textblock":
                    textblocks = find_textblocks(field[subfield], textblocks, current = current + ">" + subfield)
    return textblocks

def process_textblock(textblock):
    utext = textblock
    lines = re.sub(r"([^\n.!?\"])\n\n",r"\1.\n\n",utext)
    lines = re.sub("\n*[ ]+", " ", lines)
    lines = lines.strip()
    return lines

def get_textblocks(textblocks):
    textblocks_dict = {}
    for textblock in textblocks:
        levels = re.split(">", textblock[0])
        key = ">".join(levels[1:])
        textblocks_dict[key] = process_textblock(textblock[1])
    return textblocks_dict

def extract_textblocks(trial_id, root_dir = "E:/bio_TM/clinical_trials_data/"):
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
    return found

compounded = json.load(open("E:/bio_TM/compounded.json", "r"))

def retrieve_sentences(group = "root>burden", 
                       indicator = "(?:(?:minute|min\\.?|hour|hr\\.?|day|week|wk\\.?|month|mo\\.?|year|yr\\.?)s? (?:\\d+-|\\d+\\.)*\\d+|(?:\\d+-|\\d+\\.)*\\d+ (?:minute|min\\.?|hour|hr\\.?|day|week|wk\\.?|month|mo\\.?|year|yr\\.?)s?)", 
                       secIDs = ["brief_summary", "detailed_description"], 
                       n = 10, 
                       seed = 42,
                       index = compounded,
                       root_dir = "E:/bio_TM/clinical_trials_data/"):
    secIDs = re.compile("|".join(secIDs))
    sentIDs = list(filter(lambda x: secIDs.search(x), index.get(group, {}).get(indicator, [])))
    random.seed(seed)
    lines = []
    
    if len(sentIDs) > n:
        sentIDs = random.sample(sentIDs, n)
    for sentID in sentIDs:
        IDparts = re.split("\.", sentID)
        textblocks = {}
        textblocks = extract_textblocks(IDparts[0], root_dir = root_dir)
        line = nltk.sent_tokenize(textblocks[IDparts[1]])
        lines.append(line)

    return lines, sentIDs

def print_parse(sentence):
    doc = nlp(sentence)
    for token in doc:
        print(token.i, token.text, token.dep_, token.pos_, token.head.text, token.head.pos_,
             [child for child in token.children])

                        
def find_start(pieces, sentence):
    for piece in pieces:
        for i in range(len(sentence) - len(piece) + 1):
            j = i + len(piece)
            segment = sentence[i:j]
            ancestors = set()
            if tuple([token.text for token in segment]) == piece:
                for token in segment:
                    token_ancestors = set([token])
                    token_ancestors = token_ancestors.union(token.ancestors)
                    if ancestors:
                        ancestors = ancestors.intersection(token_ancestors)
                    else:
                        ancestors = ancestors.union(token_ancestors)
                for token in ancestors:
                    if not set(token.children).intersection(ancestors):
                        yield token, segment
                        
                        
def extract_temporal(sentence, matches):
    sentence = nlp(sentence)
    indicator_tokens = set()
    for match in matches:
        indicator_tokens.add(tuple([token.text for token in nlp(match)]))
    fragment = set()
    for start, segment in find_start(indicator_tokens, sentence):
        head = start
        for token in head.subtree:
            fragment.add(token)
        while head.pos_ != "VERB" and head.dep_ != "ROOT":
            head = head.head
            fragment.add(head)
    return sorted(fragment, key = lambda x: x.i)

examples, _ = retrieve_sentences(seed = 4)
indicator = "(?:(?:minute|min\\.?|hour|hr\\.?|day|week|wk\\.?|month|mo\\.?|year|yr\\.?)s? (?:\\d+-|\\d+\\.)*\\d+|(?:\\d+-|\\d+\\.)*\\d+ (?:minute|min\\.?|hour|hr\\.?|day|week|wk\\.?|month|mo\\.?|year|yr\\.?)s?)"
for example in examples:
     example = str(example)
     print("Sentence: ", example)
     print("Extract: ", extract_temporal(example, re.compile(indicator).findall(example)))
     print("")
