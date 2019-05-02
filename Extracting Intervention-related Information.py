# -*- coding: utf-8 -*-
"""
Created on Thu May  2 13:31:38 2019

@author: dell
"""
from __future__ import division
from collections import defaultdict, Counter
import xml.etree.ElementTree as ET
import re, json, nltk, random, spacy

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


def print_parse(sentence):
    doc = nlp(sentence)
    for token in doc:
        print(token.text, token.dep_, token.pos_, token.head.text, token.head.pos_,
             [child for child in token.children])

compounded = json.load(open("E:/bio_TM/compounded.json", "r"))     
   
def retrieve_sentences(group = "root>burden", 
                       indicator = "", 
                       secIDs = ["brief_summary", "detailed_description"], 
                       n = 10, 
                       seed = 42,
                       index = compounded,
                       root_dir = "E:/bio_TM/clinical_trials_data/"
                      ):
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

nlp = spacy.load('en_core_web_sm')

def extract_intervention(sentence, matches):
    sentence = nlp(sentence)
    indicator_tokens = set()
    for match in matches:
        indicator_tokens.add(tuple([token.text for token in nlp(match)]))
    
    for start, segment in find_start(indicator_tokens, sentence):
        verb_fragment = set()
        subj_fragments = []
        obj_fragments = []
        children = set()
        subjs = []
        objs = []
        
        head = start
        verb_fragment.add(head)
        
        for child in head.children:
            children.add(child)

        while children:
            child = children.pop()
            if re.search("subj", child.dep_):
                subjs.append(child)
            elif re.search("dobj", child.dep_):
                objs.append(child)
            elif re.search("neg|aux", child.dep_):
                verb_fragment.add(child)
            else:
                for grandchild in child.children:
                    children.add(grandchild)
                    
        verb_fragment = sorted(list(verb_fragment), key = lambda x: x.i)
        
        for subj in subjs:
            subj_fragment = set()
            for token in subj.subtree:
                if token.pos_ == "NOUN" or token.pos_ == "PROPN":
                    subj_fragment.add(token)
            parent = subj
            while parent != head:
                subj_fragment.add(parent)
                parent = parent.head
            subj_fragment = sorted(list(subj_fragment), key = lambda x: x.i)
            subj_fragments.append(subj_fragment)
        
        for obj in objs:
            obj_fragment = set()
            for token in obj.subtree:
                if token.pos_ == "NOUN" or token.pos_ == "PROPN":
                    obj_fragment.add(token)
            parent = obj
            while parent != head:
                obj_fragment.add(parent)
                parent = parent.head
            obj_fragment = sorted(list(obj_fragment), key = lambda x: x.i)
            obj_fragments.append(obj_fragment)
            
        fragment = [subj_fragments, verb_fragment, obj_fragments]
        yield fragment
        

def print_intervention_output(sentence, verb):
    print(sentence)
    print("")
    outputs = list(extract_intervention(sentence, [verb]))
    
    for output in outputs:
        subject_string = ""
        object_string = " "

        if len(output[0]) > 1:
            for subject in output[0]:
                words = [token.text for token in subject]
                single_subject_string = " ".join(words)
                subject_string = subject_string + single_subject_string + ", "
            subject_string = subject_string[:-2] + " "
        elif len(output[0]) == 1:
            subject_string += " ".join([token.text for token in output[0][0]]) + " "

        verb_string = " ".join([token.text for token in output[1]])

        if len(output[2]) > 1:
            for objet in output[2]:
                words = [token.text for token in objet]
                single_object_string = " ".join(words)
                object_string = object_string + single_object_string + ", "
            object_string = object_string[:-2]
        elif len(output[2]) == 1:
            object_string = object_string + " ".join([token.text for token in output[2][0]])

        print(subject_string + verb_string + object_string)
        print("")
    print("____________________________________________\n")
    

examples, sentIDs = retrieve_sentences(group = "root>burden>active period", indicator = u"receive", n = 20, seed = 42)
for example in examples:
    for sentence in example:
        print_intervention_output(sentence, "receive")