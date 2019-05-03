# clinical_trial
The exploration of clinical study records from ClinicalTrial.gov
### getting the data
Download all the study records data from ClinicalTrial.gov by the link: https://clinicaltrials.gov/AllPublicXML.zip.
### Extracting_Sentences_of_Interest.py
It aims to extract sentences that may contain valuable information about patients.

By load_trial() function, we can pass in a trial ID number and the location of the directory where we extracted the data, then return a Python dictionary object holding the contents of the XML file. The process_node() function can recursively navigate the XML tree and copy the data over to a dictionary.

find_textblocks() identifies the fields containing unstructured text in the loaded trial dictionary and returns a list of tuples containing the name of the field and the text in it. process_textblock() takes a single textblock and performs a little clean-up on the text data by removing extra whitespace. get_textblocks() creates a dictionary from the list of tuples with the field names as keys and the textblocks as values. And extract_textblocks() is a wrapper for all of these functions.

create_regex_format() compiles a regular expression (regex) from the defined indicator pattern. Compiling this regex allows us to match indicators at any position—beginning, middle, or end. find_indicators() processes the schema to generate a dictionary containing the compiled regex from each indicator, as well as the hierarchical information about that indicator. match_indicators() is where we use these regexes to find sentences where a match exists and put these matched fragments, along with a sentence ID number, in a dictionary.
### Extracting_Temportal_Information.py
It aims to explore the use of a syntactic parsing-based strategy to extract information that is potentially about scheduling for patients from unstructured trial description text.

Use generate_index.py to generate file named compounded. I didn't upload it because it's too large.

the function retrieve_sentences() can return a collection of n sentences matching the criteria for group, indicator, and section. The default indicator is the combination of all temporal patterns, capable of picking out mentions of minutes, hours, days, weeks, months, and years. 

print_parse() can transform a sentence into a list of spacy tokens using the spacy NLP engine. extract_temporal() can extract a readable fragment of the sentence containing the temporal pattern. find_start() initiaite the tree navigation in extract_temporal.
### Extracting Intervention-related Information.py
It aims to use specific verbs as indicators, and from a sentence's syntactical dependency tree, follow the verb to collect sentence fragments that might contain potentially useful information related to interventions.

extract_intervention() use specified verbs such as "receive" as indicators to extract intervention-related information.

print_intervention_output() prints the output in a readable format.

Thus, some intervention-related information that might be useful are collected.
