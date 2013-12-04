#!/usr/bin/env python
# :noTabs=true:

"""
remotely download from SwissModel

internal method "download_models_from_swissmodel" is the real winner
other methods are for reloading afterwards (if necessary)

...I know less about the options/details here etc., just a quick script for
downloading remotely

ugh...more messy comparative model data mining...because these people don't
make the most intuitive interfaces :(

yay! however ugly, it works! and a LOT faster than ModBase
    -do proper urllib

    -do proper taring (change filenames)
    -avoid pwd swapping etc.
    -extract more data! a lot of fun stuff on the initial swissmodel page

todo:
    rework: download, then extract!
    add in quality features...

Author: Evan H. Baugh
"""

################################################################################
# IMPORT

# common modules
import os
import re
import tarfile
import urllib2
import shutil
import optparse

# bigger modules

# custom modules
#from helper import create_directory , move_file , copy_file
#from quality import evaluate_quality , calculate_expected_rmsd_from_sequence_identity
# place these simple helper scripts here instead, stand alone!

#from secstruct import morph_secstruct
# this is just for customizing secondary structure characters
# in many cases, an extreme luxury feature...so lets abandon it for this

################################################################################
# HELPER METHODS

# not very interesting, just ignore these
# here because I like defining my own ways for these utilities which are in
# another module...but you only care about ModBase stuff right? so
# just put a copy here, that simple...

# helper for creating a directory, checks and delets existing name
def create_directory( dir_name , tagline = ' to sort the data' ):
    """
    Creates the directory  <dir_name>
    
    WARNING: this will delete the directory and its contents if it already
    exists!
    
    Optionally output something special in  <tagline>
    """
    # check if it exists
    print 'Creating a new directory ' + os.path.relpath( dir_name ) + tagline
    if os.path.isdir( dir_name ):
        print 'a directory named ' + os.path.relpath( dir_name ) + ' already exists, deleting it now...'
        shutil.rmtree( dir_name )
    os.mkdir( dir_name )

# copy helper
def copy_file( filename , destination , display = False ):
    """
    Copy  <filename>  to/into  <destination>
    
    just a cp wrapper...what?
    """
    if display:    # optional
        if os.path.isdir( destination ):
            print 'placing a copy of ' + os.path.relpath( filename ) + ' into the ' + os.path.relpath( destination ) + ' directory'
        elif os.path.isfile( destination ):
            print 'copying ' + os.path.relpath( filename ) + ' to ' + os.path.relpath( destination )
    shutil.copy( filename , destination )

# move helper
def move_file( filename , destination , display = False ):
    """
    Move  <filename>  to/into  <destination>
    
    just a mv wrapper...what?
    """
    if display:    # optional
        if os.path.isdir( destination ):
            print 'moving ' + os.path.relpath( filename ) + ' into the ' + os.path.relpath( destination ) + ' directory'
        elif os.path.isfile( destination ):
            print 'renaming ' + os.path.relpath( filename ) + ' to ' + os.path.relpath( destination )
    shutil.move( filename , destination )

################################################################################
# METHODS

# ...meh!
def download_models_from_swissmodel( query , out_directory = 'swissmodel_models' , root_filename = '' ):
    """
    REQUIRES INTERNET CONNECTION    

    Returns "details" on the models for  <query>  in SwissModel
    write results to  <out_directory>  with the base  <root_filename>

    ...
    """
    # url
    url = 'http://swissmodel.expasy.org/repository/smr.php'
    
    # format the search query
    print 'searching swissmodel for \"' + query +'\"...'
    url1 = url + '?sptr_ac=' + query
    print url1

    # soooooo
    # first we need to look at the page...so we can find the proper mapping
    # this is SO BAD but there doesn't seem to be an easy way around this...?
    raw_stream = urllib2.urlopen( url1 )
    raw_stream = raw_stream.read()
    
    # find the download link
    target = re.findall( 'href.+download model.tgz' , raw_stream )

    # exit condition, no link found    
    if not target:
        # exit condition
        print 'no models exist in Swissmodel for this protein...'
        return None
    else:
        target = target[0]    # should only be 1, the first hit
    target = target[6:-20]    #...oh I hope these are stable...
    print 'found where the real models are stored, downloading now...'
    
    # grab the actual download
    url2 = url + target
    raw_stream = urllib2.urlopen( url2 )
    raw_stream = raw_stream.read()
    
    # optionally make a new directory
    if out_directory:
        create_directory( out_directory , ' to store the models, alignments, etc.' )

        # so bad...
        # for safe later writing...
        out_directory += '/'

    # defaults for writing files
    if not root_filename:
        root_filename = 'swissmodel_' + query

    # write the loaded stuff to a file
    filename = out_directory + root_filename + '.tgz'
    print 'writing the compressed download to ' + filename
    f = open( filename , 'w' )
    f.write( raw_stream )
    f.close()

    # decompress/load the contents of the tar file
    print 'extracting individual files...'
    results = tarfile.open( filename , 'r' )
    
    # find the annoying database str...
    files = results.getmembers()
    database_str = [i.path[:-10] for i in files if i.path[-10:] == '.templates']
    if database_str:
        database_str = database_str[0]
        # assume only 1 will be here and it will match the ".templates" file...
    
    # extract the results
    for i in results:
        results.extract( i , path = out_directory )
    
    # rename the files
    new_files = []
    for i in files:
        new_filename = i.path
        if database_str and database_str in i.path:
            new_filename = out_directory + i.path.replace( database_str , root_filename )
            move_file( out_directory + i.path , new_filename )

        new_files.append( new_filename )

    # if no database_str found, learn it from the files
    # an exception case thats been causing problems
    if not database_str:
        print 'Odd record!!! FYI this is one of those weird exceptions...'
        # consider the filenames
        database_str = [i.split( '_' )[0] for i in os.listdir( out_directory ) if 'MODEL.txt' in i][0]
        # assume only 1 ID esists
        for i in os.listdir( out_directory ):
            if database_str in i:
                new_filename = out_directory + i.replace( database_str , root_filename )
                move_file( out_directory + i , new_filename )

        # update "new_files"
        new_files = [out_directory + i.replace( database_str , root_filename ) for i in new_files]

    # determine the number of models
    model_summaries = [i for i in new_files if i[-10:] == '.MODEL.txt']
    model_details = {}
    for i in model_summaries:
        # get the model ID/number...
        model_number = i.split( '_' )[-1]    # the last one...
        model_number = int( model_number[:model_number.find( '.' )] )
        
        # extract the summary
        model_details[model_number] = extract_model_details_from_swissmodel_summary( i )
        model_details[model_number]['coordinates'] = os.path.abspath( i[:-10] + '.pdb' )

    # convert to a list
    temp = model_details.keys()
    temp.sort()
    model_details = [model_details[i] for i in temp]

    # find the "best" model
    print '\nevaluating the \"best\" model by comparing:\n\t1. sequence identity\n\t2. target length'
    best_score = max( [i['sequence identity'] for i in model_details] )
    matches = [i for i in model_details if i['sequence identity'] == best_score]
    if len( matches ) > 1 and sum( [not i['target length'] == matches[0]['target length'] for i in matches[1:]] ):
        best_score = max( [i['target length'] for i in model_details] )
        matches = [i for i in model_details if i['target length'] == best_score]
        # add "best" to model details?
   
    # debug output
    if len( matches ) > 1:
        print 'multiple models are \"equally the best\":'
        for i in matches:
            print '\t'+ i['coordinates']
        print 'copying the first one to best_model.pdb'
    else:
        print 'best model: ' + matches[0]['coordinates']
    # move it to a indicative filename
    copy_file( matches[0]['coordinates'] , out_directory + '/best_model.pdb' )

    return model_details

############################
# extract additional details

# returns a summary of a model by scanning the MODEL.txt files
def extract_model_details_from_swissmodel_summary( swissmodel_summary_filename ):
    """
    Returns a dict of summary lines in  <swissmodel_summary_filename>
    Thank you SwissModel, a straight-forward file format! 1 empty header, then
    a bunch of "=" separated terms
    filter down to the useful ones...
    """
    # load the lines
    f = open( swissmodel_summary_filename , 'r' )
    lines = [i.strip().split( '=' ) for i in f.xreadlines() if i.strip()]
    f.close()

    # I despise doing this...but here, exposure of this is MORE confusing
    # than obfuscating/hiding it in the code
    # similar to other parsing jobs, this set is unlikely to change
    # and pretty much useless as "quicker" modular changing
    acceptable_terms = {
        'sequence' : 'sequence' ,
        'align_identity' : 'sequence identity' ,
        'target_start' : 'target start' ,
        'target_stop' : 'target end' ,
        'template_start' : 'template start' ,
        'template_stop' : 'template end' ,
        'align_evalue' : 'alignment evalue' ,
        'model_start' : 'model start' ,    # what do these mean?
        'model_end' : 'model end' ,        # what do these mean?
        'template' : 'template' ,
        'chain' : 'template chain' ,
        'modeling_method' : 'program' ,
        'model_date' : 'model date' ,
        'align_method' : 'alignment program'
        }
    
    # feeling lazy, simple filter
    model_details = { 'model date' : '???' }
    for i in lines:
        key = i[0].strip()
        if key in acceptable_terms.keys():
            model_details[acceptable_terms[key]] = i[1].strip()

    # post processing? make it like modbase?
    if 'target start' in model_details.keys() and 'target end' in model_details.keys():
        model_details['target coverage'] = [int( model_details['target start'] ) , int( model_details['target end'] )]
        model_details['target length'] = int( model_details['target end'] ) - int( model_details['target start'] ) + 1
        # and clean up...
        del model_details['target start']
        del model_details['target end']
    if 'template start' in model_details.keys() and 'template end' in model_details.keys():
        model_details['template coverage'] = [int( model_details['template start'] ) , int( model_details['template end'] )]
        model_details['template length'] = int( model_details['template end'] ) - int( model_details['template start'] ) + 1
        # and clean up...
        del model_details['template start']
        del model_details['template end']
    # what is the difference between the model and the target!? what?
    if 'model start' in model_details.keys() and 'model end' in model_details.keys():
        model_details['model coverage'] = [int( model_details['model start'] ) , int( model_details['model end'] )]
        model_details['model length'] = int( model_details['model end'] ) - int( model_details['model start'] ) + 1
        # and clean up...
        del model_details['model start']
        del model_details['model end']
    
    if 'template' in model_details.keys():
        model_details['template'] = model_details['template'].upper()
    if 'alignment evalue' in model_details.keys():
        model_details['alignment evalue'] = float( model_details['alignment evalue'] )
    if 'sequence identity' in model_details.keys():
        model_details['sequence identity'] = float( model_details['sequence identity'] )
    
    # get covered sequence
    if 'sequence' in model_details.keys() and 'target coverage' in model_details.keys():
        model_details['target sequence'] = model_details['sequence'][model_details['target coverage'][0] - 1:model_details['target coverage'][1]]
    
    return model_details

# ...um...useful resource for finding ss predictions? etc.?
def extract_modeling_details_from_swissmodel_hhm( swissmodel_hhm_filename ):
    """
    Returns the secstruct string from HHPred recorded in the SwissModel HHM
    file, a part of the modeling record
    
    ...currently this is the only "useful" information in this file, also
    contains alignment information, but this appears to be a component of the
    template search...so not too useful in analyses
    could extract confidence of the secstruct prediction...or consensus...naw...
    
    currently uses HHPred ss character str...should convert...
    """
    # load the file
    f = open( swissmodel_hhm_filename , 'r' )
    in_secstruct_str = False
    secstruct_str = ''
    for line in f.xreadlines():
        # scanning for the start
        if line[:8] == '>ss_pred':
            # start reading the secstruct str
            in_secstruct_str = True
            continue    # don't add it!
        elif line[:8] == '>ss_conf':
            # exit, no reason to keep reading...
            break

        if in_secstruct_str:
            secstruct_str += line.strip()
    f.close()
    
    # convert the str?
    # oh well...
#    secstruct_str = morph_secstruct( secstruct_str )
    # not supported in this version of the script
    
    return secstruct_str

# returns a summary of additional details found on the swissmodel page
def extract_additional_details_from_swissmodel_html_page( swissmodel_html_str ):
    None
    # model quality
    # alignments/summary...already in downloads
    # oligomeric, ligands and complex/interactions etc.
    
    # ...not really useful?...mane, this is tedious to automate and
    # the "quality" is completely ???

################################################################################
# MAIN

if __name__ == '__main__':
    # parser object for managing input options
    parser = optparse.OptionParser()
    # essential data
    parser.add_option( '-q' , dest = 'query' ,
        default = '' ,
        help = 'the (UniProt) ID to download from SwissModel' )

    parser.add_option( '-o' , dest = 'out_directory' ,
        default = 'swissmodel_models' ,
        help = 'the name of the directory to create, for storing the downloads' )
    parser.add_option( '-r' , dest = 'root_filename' ,
        default = '' ,
        help = 'the \"root\" for naming the files, defaults to the query' )

    (options,args) = parser.parse_args()

    # check inputs
    # no edits/modifications
    # kinda silly, but I do this as "my style", easy to modify cleanly
    query = options.query
    out_directory = options.out_directory
    root_filename = options.root_filename

    # choose the default method to run
    download_models_from_swissmodel( query , out_directory = out_directory ,
        root_filename = root_filename )


