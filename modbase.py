#!/usr/bin/env python
# :noTabs=true:

"""
remotely download from modbase

so...some quirks, there seem to be multiple "user levels" including a
distinction between "public" and "academic"...this seems odd to me so for now
just looking at the public interface...or academic...as it seems to be?

reference the modeller.py script for scoring info...sorta...

there may be other data that is extractable from ModBase...however this is a bit
obfuscated, for now, just focus on downloading the models and alignments

rework: download, then extract!
add in quality features...

Author: Evan H. Baugh
"""

################################################################################
# IMPORT

# common modules
import os
import urllib2
from xml.dom.minidom import parse as xml_parse

# bigger modules

# custom modules
from helper import create_directory , copy_file

################################################################################
# METHODS

# woohoo!
def download_models_from_modbase( query ,
        out_directory = 'modbase_models' , root_filename = '' ,
        dataset = '' , get_alignment = True , write_summary = True ,
        display = True ):
    """
    REQUIRES INTERNET CONNECTION    

    Returns "details" on the models for  <query>  in ModBase
    write results to  <out_directory>  with the base  <root_filename>
    
    Optionally request models from a specific  <dataset>
    Optionally  <get_alingment>  too (as PIR file)
    Optionally  <display>  a summary of the results
    Optionally  <write_summary>  of the models (human readable, also displays)
    
    ModBase documentation claims that the interface can accept:
        databaseID      database ID, let's use UniProt
        dataset         a particular ModBase run?
        modelID         same?
        seqID           same?

        dataset         the ModWeb JobID...
        type            "model" or "alignment", this method handles this
    and that any of the first 4 is enough to identify the target (?)
    ...for simplicity, let's just look using UniProt IDs as "databaseIDs"
    
    apparently to use "non-public" access additional work must be done
    (something about a "cookies.txt" file, though this seems specific to "wget",
    may be able to pass in user/password as "modbase_user" and "modbase_passwd")
    
    uses xml.dom.minidom to parse the HTML returned...this may not be kosher...
    but it works...and is easier than using htmllib or sgmllib...(?)
    """
    # url
    url = 'http://salilab.org/modbase/retrieve/modbase'
    
    # format the search query
    print 'searching modbase for \"' + query +'\"'
    url += '?databaseID=' + query
    # currently unused...so why put it here?
    #for i in search_options.keys():
    #    url += '&' + i +'='+ search_options[i]
    
    # the dataset
#    if not 'dataset' in search_options.keys() and dataset:
    if dataset:
        url += '&dataset=' + dataset

    # go get the results
    print 'obtaining model results from:\n\t' + url
    raw_stream = urllib2.urlopen( url  + '&type=model' )    
    print 'finished downloading models, summarizing the results...'
        
    # parse the results
    results = xml_parse( raw_stream )

    # check if empty
    if not len( results.toxml() ) > 100:    # ahhh! I hate arbitrary numbers!!!
        print 'no models exist in ModBase for this protein...'
        return {}
    
    # get the ids
    #ids = get_str_from_xml_tag( results , 'model_id' )
    # no need, in the header of the model
    
    # get the models
    models = get_str_from_xml_tag( results , 'content' )
        
    # extract the details
    details , text = get_modbase_model_details( models , display or write_summary , export = True )
        
    # defaults for writing files
    if not root_filename:
        root_filename = 'modbase_' + query
    
    # optionally write the models
    if out_directory:
        create_directory( out_directory , ' to store the models as PDB files' )
        print 'writing the downloaded models to ' + out_directory
        count = 1
        filenames = []
        for i in models:
            # write it
            filename = out_directory + '/' + root_filename + '_model_' + str( count ) + '.pdb'
            filenames.append( os.path.abspath( filename ) )

            # write the alignment
            f = open( filename , 'w' )
            f.write( i.strip() )
            f.close()
            count += 1
        
        # change this in this case
        models = filenames
        
        # SOOO HACKY!!!!
        # for later safety...
        out_directory += '/'

    # optionally grab the alignment too
    if get_alignment:
        print 'also downloading the alignments...'
        raw_aln_stream = urllib2.urlopen( url  + '&type=alignment' )

        # parse the results
        aln_results = xml_parse( raw_aln_stream )
        
        # get the files
        aln_results = aln_results.getElementsByTagName( 'alignmentfile' )
        
        # ...for now, just get the text itself
        # don't worry about the other details in the XML file
        print 'writing the alignments as PIR files...'
        count = 1
        for i in aln_results:
            i = get_str_from_xml_tag( i , 'content' )[0]    # just 1, always the first
            
            # if out_directory is empty...this will just do as we want
            filename = out_directory + root_filename + '_model_' + str( count ) + '_alignment.pir'
            f = open( filename , 'w' )
            f.write( i )
            f.close()
            
            # convert them?
            # doesn't seem to load these "pir" files...? :(
            
            # save in the details?
            details[count - 1]['alignment'] = i
            
            count += 1
    
    # put the models (filenames) into details...cleaner output, just 1 dict
    for i in xrange( len( models ) ):
        details[i]['coordinates'] = models[i]
    
    # find the "best" model
    temp = '\nevaluating the \"best\" model by comparing:\n\t1. sequence identity\n\t2. model score\n\t3. target length'
    print temp
    text += temp +'\n'
    best_score = max( [i['sequence identity'] for i in details] )
    matches = [i for i in details if i['sequence identity'] == best_score]
    if len( matches ) > 1 and sum( [not i['model score'] == matches[0]['model score'] for i in matches[1:]] ):
        # find the best model score
        best_score = max( [i['model score'] for i in details] )
        matches = [i for i in details if i['model score'] == best_score]
        
        if len( matches ) > 1 and sum( [not i['target length'] == matches[0]['target length'] for i in matches[1:]] ):
            best_score = max( [i['target length'] for i in details] )
            matches = [i for i in details if i['target length'] == best_score]
   
    # debug output
    if len( matches ) > 1:
        temp = 'multiple models are \"equally the best\":'
        print temp
        text += temp +'\n'
        for i in matches:
            temp = '\t'+ i['coordinates']
            print temp
            text += temp +'\n'
        temp = 'copying the first on to best_model.pdb'
        print temp
        text += temp +'\n'
    else:
        temp = 'best model: ' + matches[0]['coordinates']
        print temp
        text += temp
    # move it to a indicative filename
    copy_file( matches[0]['coordinates'] , out_directory + '/best_model.pdb' )

    # optionally write a summary file
    if write_summary:
        # if out_directory is empty...this will just do as we want
        filename = out_directory + root_filename + '_summary.txt'
        f = open( filename , 'w' )
        f.write( text )
        f.close()
    
    # just the details, has everything else...
    return details

# very hacky wrapper
def get_str_from_xml_tag( xml_object , tag ):
    """
    So...I don't have time to learn proper XML parsing with the Python "xml"
    library right now and this approach works...so yeah
    
    simply return a list of str for the target  <tag>  in  <xml_object>
    """
    # get it
    results = xml_object.getElementsByTagName( tag )
    
    # convert to string
    L = len( tag )
    results = [i.toxml()[L + 2:-(L + 3)].strip() for i in results]
    
    return results

# useful simple text parsing
def extract_model_details_from_modbase_header( modbase_model_text ):
    """
    Returns a dict of the model details from  <modbase_model_text>
    this includes the PDB template, coverage details (always continuous),
    and alignment/modeling details
    """
    # setup defaults, cleaner display
    details = {
        'model' : '?' ,
        'organism' : '?' ,
        'experiment' : '?' ,
        'method' : '?' ,
        'program' : '?' ,

        'sequence identity' : 0 ,
        'model score' : 0 ,
        'evalue' : 0 ,

        'template' : '?' ,
        'template chain' : '?' ,
        'template coverage' : [] ,
        'target length record' : 0 ,
        'target coverage' : [] ,
        'template length' : 0 ,
        'target length' : 0 ,
        
        'ModPipe run' : '?' ,
        'modelID' : '?' ,
        'alignmentID' : '?'
        }
        
    # over the lines
    for i in modbase_model_text.split( '\n' ):
        if i[:4] == 'ATOM':
            # done! end of the header
            break
        elif i[:6] == 'HEADER':
            details['model'] = str( i.replace( 'HEADER' , '' ).strip() )
        #elif i[:5] == 'TITLE':    # ...uh, the ones I looked at, this was useless...
        elif i[:6] == 'SOURCE':
            details['organism'] = str( i.replace( 'SOURCE' , '' ).strip() )
        #elif i[:6] == 'AUTHOR':    # don't care about authors for now...
        elif i[:10] == 'REMARK 220':
            i = str( i.replace( 'REMARK 220' , '' ).strip() )
            
            # keep sorting...
            if i[:16] == 'EXPERIMENT TYPE:':
                details['experiment'] = str( i.replace( 'EXPERIMENT TYPE:' , '' ).strip() ).capitalize()
            elif i[:7] == 'METHOD:':
                details['method'] = str( i.replace( 'METHOD:' , '' ).strip() ).capitalize()
            elif i[:8] == 'PROGRAM:':
                details['program'] = str( i.replace( 'PROGRAM:' , '' ).strip() )

            elif i[:18] == 'SEQUENCE IDENTITY:':
                # as fraction please
                details['sequence identity'] = float( i.replace( 'SEQUENCE IDENTITY:' , '' ).strip() )/100
            elif i[:12] == 'MODEL SCORE:':
                # as float
                details['model score'] = float( i.replace( 'MODEL SCORE:' , '' ).strip() )
            elif i[:7] == 'EVALUE:':
                # as float
                details['evalue'] = float( i.replace( 'EVALUE:' , '' ).strip() )

            elif i[:13] == 'TEMPLATE PDB:':
                details['template'] = str( i.replace( 'TEMPLATE PDB:' , '' ).strip().upper() )
            elif i[:15] == 'TEMPLATE CHAIN:':
                details['template chain'] = str( i.replace( 'TEMPLATE CHAIN:' , '' ).strip() )
            elif i[:15] == 'TEMPLATE BEGIN:':
                details['template coverage'].append( int( i.replace( 'TEMPLATE BEGIN:' , '' ).strip() ) )
            elif i[:13] == 'TEMPLATE END:':
                details['template coverage'].append( int( i.replace( 'TEMPLATE END:' , '' ).strip() ) )

            elif i[:14] == 'TARGET LENGTH:':
                details['target length record'] = int( i.replace( 'TARGET LENGTH:' , '' ).strip() )
            elif i[:13] == 'TARGET BEGIN:':
                details['target coverage'].append( int( i.replace( 'TARGET BEGIN:' , '' ).strip() ) )
            elif i[:11] == 'TARGET END:':
                details['target coverage'].append( int( i.replace( 'TARGET END:' , '' ).strip() ) )

            elif i[:12] == 'MODPIPE RUN:':
                details['ModPipe run'] = str( i.replace( 'MODPIPE RUN:' , '' ).strip() )
            elif i[:17] == 'MODPIPE MODEL ID:':
                details['modelID'] = str( i.replace( 'MODPIPE MODEL ID:' , '' ).strip() )
            elif i[:21] == 'MODPIPE ALIGNMENT ID:':
                details['alignmentID'] = str( i.replace( 'MODPIPE ALIGNMENT ID:' , '' ).strip() )

    # for own sanity
    details['template coverage'].sort()
    details['template length'] = details['template coverage'][1] - details['template coverage'][0] + 1
    details['target coverage'].sort()
    details['target length'] = details['target coverage'][1] - details['target coverage'][0] + 1

    return details

# silly interactive method
def display_modbase_model_details( details , include_run_details = False , export = False ):
    """
    Displays a summary of the ModBase model  <details>
    
    Optionally  <include_run_details>
    Optionally  <export>  the summary text
    """
    # check the input
    if isinstance( details , str ):
        # assume it just needs to be parsed out
        details = extract_model_details_from_modbase_header( details )
    
    # exit condition
    if 'FAIL' in details.keys():
        text = details['model'] +'\n'
        text += 'FAILED: ' + details['FAIL'] +'\n'
        
        print text[:-1]
        
        return text[:-1]
    
    # make the string
    text = details['model'] +'\n'
    text += 'covering: ' + str( details['target coverage'][0] ) +'-'+ str( details['target coverage'][1] ) +' ('+ str( details['target length'] ) + ' positions)\n'
    
    text += '\nfrom: ' + details['template'] +' '+ details['template chain'] +' (' + str( details['template coverage'][0] ) +'-'+ str( details['template coverage'][1] ) +') from ' + details['organism'] +'\n'
    
    text += '\nsequence identity: ' + str( details['sequence identity'] )[:6] + ' (evalue ' + str( details['evalue'] )[:6] +')\n'
    text += 'model score: ' + str( details['model score'] )[:6] +'\n'
    
    # optionally include the run details
    if include_run_details:
        text += '\n'
        text += details['experiment'] + ' by ' + details['method'] + ' using ' + details['program'] +'\n'
        text += 'from: ' + details['ModPipe run'] +'\n'
    
    print text[:-1]
    
    # optionally return the text
    if export:
        return text[:-1]

# simple wrapper
def get_modbase_model_details( models , add_model_numbers = True , display = True , export = False ):
    """
    Returns the details of the model text  <models>
    
    Optionally  <display>  the details
    Optionally  <export>  the details AND a summary str
    
    ...this is quite messy, but if the interpreter closes, the data must be
    extracted from the downloaded files, hence replicate the naming conventions
    used...
    """
    # get the details
    # add hacky error check here...
    permissive_models = []
    for i in models:
        try:
            model = extract_model_details_from_modbase_header( i )
            permissive_models.append( model )
        except:
            print 'HEY! very rare bug found! ModBase model based on icode regions...not today friends...'
            permissive_models.append( {
                'FAIL' : 'model indices contain and icode...this makes re-aligning VERY complicated...' ,
                'model' : '???' , 
                'sequence identity' : 0 ,
                'model score' : 0 ,
                'target length' : 0 ,    # more hacking...ugh
                'target coverage' : [0 , 0]
                } )
            #continue
    details = permissive_models
    
    # optionally add the model # into the name
    if add_model_numbers:
        for i in xrange( len( details ) ):
            details[i]['model'] = '#' + str( i + 1 ) +', '+ details[i]['model']

    # make text for display/writing
    # optionally display the model summay...because of how this is setup...
    if display:
        text = ''
        for i in xrange( len( details ) ):
            temp = '='*80
            print temp
            text += temp +'\n'
            text += display_modbase_model_details( details[i] , export = True ) +'\n'
        print    # clearner...

    # optionally return the details and a text summary...hacky...
    if export:
        return details , text
    return details

# messy...based on default naming...
def determine_modbase_models_from_modbase_directory( query , out_directory = 'modbase_models' , root_filename = '' ):
    """
    Returns the "expected" model filenames from  <out_directory>  downloaded
    from ModBase
    
    ...this is quite messy, but if the interpreter closes, the data must be
    extracted from the downloaded files, hence replicate the naming conventions
    used...
    """
    # defaults for written files
    if not root_filename:
        root_filename = 'modbase_' + query
    if not out_directory:
        out_directory = './'    # here!
    
    # ta da!
    return [i for i in os.listdir( out_directory ) if root_filename + '_model_' in i and i[-4:] == '.pdb']

# ugh...
def extract_modbase_model_details_from_modbase_directory( query , out_directory = 'modbase_models' , root_filename = '' , display = True ):
    """
    Returns the model details of ModBase models for  <query>  in
    <out_directory>
    
    Optionally specify the  <root_filename>  from when the files were written
    Optionally  <display>  the details (for interactive use)
    
    ...this is quite messy, but if the interpreter closes, the data must be
    extracted from the downloaded files, hence replicate the naming conventions
    used...
    """
    # get the models
    models = determine_modbase_models_from_modbase_directory( query , out_directory , root_filename )

    # :( hacky...but must be sure
    # rather than an if...
    out_directory = out_directory.strip( '/' ) +'/'

    # load the text
    model_text = ['']*len( models )
    for i in xrange( len( models ) ):
        f = open( out_directory + models[i] , 'r' )
        model_text[i] = f.read()
        f.close()
    
    # extract the details
    details = get_modbase_model_details( model_text , display )
    
    # put the models into details...cleaner output, just 1 dict
    for i in xrange( len( models ) ):
        details[i]['coordinates'] = models[i]
        details[i]['model'] = '#' + str( i + 1 ) +', '+ details[i]['model']
        
        aln_filename = out_directory + models[i][:-4] + '_alignment.pir'
        if os.path.exists( aln_filename ):
            f = open( aln_filename , 'r' )
            details[i]['alignment'] = f.read()
            f.close()
    
    return details

################################################################################
# MAIN

if __name__ == '__main__':
    None


