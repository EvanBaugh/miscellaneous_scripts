#!/usr/bin/env python
# :noTabs=true:

"""
Methods for cleaning and parsing PDB files

Most importantly, the process_pdb method does a lot to clean PDB files
from RCSB

Requires:
    
    Biopython

Uses:

    paths.py
    settings.py
    biopython_settings.py
    seq_basics.py
    
Author: Evan H. Baugh
"""

################################################################################
# IMPORT

# common modules
import optparse    # for commandline
import os
import shutil

# bigger modules
from Bio.PDB import PDBIO
from Bio.PDB import PDBParser
from Bio.PDB import PPBuilder    # no longer used, much faster way to do this
#from Bio.PDB import Select    # no longer used...kinda hard to use
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# custom modules
from helper import get_root_filename , create_directory , copy_file
from settings import SEQFORMAT , SEQFORMAT_EXTENSION_MAP , NUCLEIC_SEQUENCE_LETTERS_MAP , NA_CODES , three2one , WATER_CONVERSION , one2three , three2three , NA_CONVERSIONS_ROSETTA
from biopython_settings import DNAAlphabet , ProteinAlphabet

from seq_basics import write_sequence , get_sequence

################################################################################
# FULL RAW PROCESSING

# 1R69 - single model, single chain
# 1A17 - another random choice for testing
# 1BUW
# 1C17
# 1JYX
# 1M2V
# 1TF6
# 2C35
# 3G3O
# 1YY8 - AB and CD, single model
# 1NMR - multiple models
# 1LR1 - multiple models AND chains
# 1VTL - protein and DNA, single model
# 1UN6 - protein and RNA, single model

# the big boy...
def process_pdb( pdb_filename , seqformat = SEQFORMAT , seqformat_extension_map = SEQFORMAT_EXTENSION_MAP , conversion = three2three , na_conversion = NA_CONVERSIONS_ROSETTA , na_alphabet = DNAAlphabet , protein_alphabet = ProteinAlphabet ):
    """
    Create a directory from  <pdb_filename>  containing relevant information
    stored in the PDB file
    
    This method behaves slightly differently for PDB files with multiple models,
    nucleic acids, duplicate complexes, etc.
    so if you are interested in the specifics, please read the source code
    
    In short, it tries to write:
        header.txt          a text file of the header lines
        numbering_map.txt   a text file showing 1-indexed PDB numbering
        clean.pdb           only ATOM lines
        hetatm.pdb          only HETATM lines, may be split by resName
        .fa                 sequences of all peptides and nucleic acids
        subdirectories      for each protein model/subunit (similar info)
    
    does not write a text file for the "trailer" (lines after the coordinates)
    
    converts lines (ATOM or HETATM) that can be converted based on  <conversion>
    (generally) and  <na_conversion>  (specific for nucleic acids, relevant
    because RNA and DNA may require different treatment...)
    !!!WARNING!!! defaults:
        CSE     CYS     converts SelenoCysteinE to Cysteine
        HYP     PRO     converts HYdroxylProline to Proline
        CYD     CYS     does NOT convert "CYsteine Disulfides to Cysteine"
        HIP     HIS     converts "HIP" to Histidine (~double protonation)
        HID     HIS     converts "HID" to Histidine (~single delta N proton)
        HIE     HIS     converts "HIE" to Histidine (~single epsilon N proton)

    todo:
    ensure hetatm conversions step illegal atoms!!!!
    alternate conformations
    convert DNA to Rosetta DNA
    convert ligands to params
    convert water to TP3 (or TP5)
    """
    # process input, optionally a list
    if isinstance( pdb_filename , list ):
        print 'Multiple PDB codes detected, processing them individually...'
        # use this list comprehension, get them all!
        filenames = [process_pdb( i , seqformat , seqformat_extension_map , conversion , na_conversion , na_alphabet , protein_alphabet ) for i in pdb_filename]
        print 'Finished the whole list, enjoy!'
        return filenames
    
    ####################
    # NEW DIRECTORY ETC.
    # get root name
    pdb_filename = os.path.abspath( pdb_filename )
    root_name = get_root_filename( pdb_filename )
    best_guess = pdb_filename
            
    # make a new directory, a whole lot is gonna go here...
    create_directory( root_name , ' to sort the data' )
    # move the pdb here
    copy_file( pdb_filename , root_name )

    # oh, and go there too
    original_dir = os.getcwd()
    os.chdir( root_name )

    # "update" the target
    pdb_filename = root_name + '/' + os.path.split( pdb_filename )[-1]
    root_name = get_root_filename( pdb_filename )
    
    ##############
    # PRE CLEANING
    # does not need to know if nucleics or not
    
    # convertions!
    # ...bad...overwrite the file!...but no filename management
    convert_pdb_resnames_to_ATOM_lines( pdb_filename , pdb_filename , root_name +'_conversion_report.txt' , conversion )
    
    # produce a PDB with just the protein lines
    best_guess = clean_ATOM_lines_from_pdb( pdb_filename )
    
    # extract numbering
    # don't bother storing the map
    extract_numbering_map_from_pdb( pdb_filename , 'numbering_map.txt' )
    
    # extract HETATM lines
    clean_HETATM_lines_from_pdb( pdb_filename )
    
    # write out alternate conformations for the cleaned file
    alternate_conformations = clean_alternate_conformations_from_pdb( best_guess )
    
    ##########################
    # HEADER PARSING
    # extract info from header

    # this information is accessible from the PDBParser header...sorta...    
    # get the number of models
    models = extract_number_of_models_from_pdb_header( pdb_filename )
    
    # get the subunit complexes
    complexes = extract_duplicate_chains_from_pdb_header( pdb_filename )
        
    # write the header (?)
    # get the header
    header = extract_header_from_pdb( pdb_filename )

    ###################
    # HUNT DOWN HETATMS
    
    # use the map in the header and extracted chemical formulas to search pubchem
    
    # get map
    
    # per hetatm type
    # get formula
    # get number of residues -> needed to interpret formula...
    # search pubchem, download best sdf if exact match and at least < atoms
    # create directory for these params etc.

    ##########################
    # ASSESS NUCLEIC SITUATION

    # HERE!
    
    # choose your fate!, removes nucleic lines
    has_nucleic = clean_nucleic_acid_lines_from_pdb( pdb_filename )
    
    # get proteins if nucleics
    if has_nucleic:
        # get a PDB of protein only, use this from now on
        print 'Scanners indicate there are nucleic acid lines in ' + os.path.relpath( pdb_filename ) + '\nSadly, a lot of toys do not play well with these so a few extra steps are required...'

        # write nucleic sequences
        temp , nucleic_types = extract_nucleic_acid_sequences_from_pdb( root_name + '.nucleic.pdb' , seqformat = seqformat , alphabet = na_alphabet , seqformat_extension_map = seqformat_extension_map )
        # care not for the sequences
        
        # make a Rosetta ready nucleic PDB!!!
        # SO BAD! overwrite!
        # BAH!!!
        na_chains = split_pdb_into_chains( root_name + '.nucleic.pdb' , 0 , True )    # just 0 model...
        for i in na_chains.keys():
            # BETTER BE IN BOTH!!!
            convert_pdb_resnames_to_ATOM_lines( na_chains[i] , na_chains[i] , 'nucleic_chain_'+ i +'_conversion_report.txt' , na_conversion[nucleic_types[i]] )
        
        # check for protein :)
        has_protein = clean_protein_lines_from_pdb( pdb_filename )
        if not has_protein:
            print 'The additional features are only available for proteins\nScanner indicate that this PDB has ONLY nucleic acids (no proteins) :(\nthe remaining methods rely on the Biopython PDBParser...and things get messy with nucleic acids\nEven so, the only feature you\' missing out on is splitting into subdirectories for each chain, and since the PDB is just nucleic acid, that isn\'t as helpful'

            # premature exit
            os.chdir( original_dir )
            return best_guess

        # change the name of the best guess to .protein.pdb
        best_guess = root_name + '.protein.pdb'
        pdb_filename = root_name + '.protein.pdb'
        
        # get the nucleic chains
        nucleic_chains = extract_chains_from_pdb( root_name + '.nucleic.pdb' )
    
    ############
    # PDB PARSER
    
    # does NOT loop over ANY nucleic acid chains!
    
    # prepare to load...
    parser = PDBParser( PERMISSIVE = 1 )
    writer = PDBIO()
    
    struct = parser.get_structure( root_name , pdb_filename )

    # verify models and chains
    temp = len( struct.child_list )    # number of models
    if not temp == models:
        print 'Huh? the PDB file header claims there are ' + str( models ) + ' models but the PDB file has ' + str( temp ) + ' models...\nUsing the ACTUAL number of models (' + str( temp ) + ')'
        models = temp

    # check from reading the CHAIN
    if not complexes:
        print 'No chain/subunit information found in the header (or no header),\nassuming all individual sequences are unique i.e. if AB and copy CD, will make A, B, C, and D instead of AB and CD'
#        complexes = temp    # unecessary, automatically happens below...

    # add all new ids
    temp = struct[0].child_dict.keys()    # it better have at least 1 model...

    
    # for the nucleic case...
    if has_nucleic:
        # HERE!
        # remove nucleic lines...
        for i in xrange( len( complexes ) ):
            for j in nucleic_chains:
                if j in complexes[i]:
                    complexes[i] = complexes[i].replace( j ,'' )
        
        # sanity check...
        complexes = [i for i in complexes if i]
        
        # assume all models contain all chains...idk how this would ever NOT occur...
        # this also produces a directory for EACH chain as the default behavior!!!
        complexes += [i for i in temp if i and not i in complexes and not i in nucleic_chains]

    else:
        # normal protein stuff
        complexes += [i for i in temp if i and not i in complexes]
    
    
    # okay...this should be figured out...but isn't that big of a deal
    # found with 1JGO
#    print complexes
#    print complexes
#    complexes = [i for i in complexes if i]
    
#    input('dd')
    
    ################################
    # CREATE AND FILL SUBDIRECTORIES
    
    # again, this step is skipped for pure nucleic acid...
                
    # exit condition, only 1 model and 1 chain
    if models > 1 or len( complexes ) > 1:
        # over the models
        for model in struct.child_dict.keys():
            # over the chains
            for complx in complexes:
#                print '='*60 + complx
                # remove nucleic subunits
                # HERE!
                if has_nucleic:
                    for chain in nucleic_chains:
                        complx = complx.replace( chain , '' )    # delete the chain from the complex
            
                # check that all members are present
                chains = struct[model].child_dict.keys()
                missing = [l for l in complx if not l in chains]
                # report this!
                if missing:
                    # add models bool for str here?
                    print 'Expected model ' + str( model + 1 ) + ' to have chains ' + complx + ' but the its missing chains ' + ', '.join( missing ) + '!'
                
                # create the new directory
                # only number if more than 1 model
                dir_name = complx + str( model + 1 )*bool( models - 1 )
                new_dir = os.path.split( root_name )[0] + '/' + dir_name
                print 'Creating the subdirectory ' + os.path.relpath( new_dir )
                os.mkdir( new_dir )

                # create a copy of the complex, only the chains of interest
                # make an empty structure
                temp = Structure( 'temp' )
                temp_model = Model( model )    # and an empty model
                temp.add( temp_model )
                
                # add the complex
                for chain in complx:
                    temp[model].add( struct[model][chain] )
                    
                    # get the chain sequence
                    seqid = dir_name + ('_model_' + str( model + 1 ))*bool( models - 1 ) + '_chain_' + chain
                    seq_filename = new_dir + '/' + os.path.split( root_name )[-1] + ('_model_' + str( model + 1 ))*bool( models - 1 ) + '_chain_' + chain + '.' + seqformat_extension_map[seqformat]
                    description = '(from model ' + str( model + 1 ) + ')'
                    temp_seq = extract_protein_sequence_from_pdb( temp , True ,    # MUST insert disorder...
                        seq_filename , seqid , description , model , chain ,
                        True , seqformat , protein_alphabet , seqformat_extension_map )

                    # also, make sure at least one copy (from the first model) is in the main dir
                    seq_filename = root_name + '_chain_' + chain + '.' + seqformat_extension_map[seqformat]
                    if not os.path.exists( seq_filename ):
                        print 'Putting a copy of the sequence in the new directory'
                        # assumes all the models have the same sequence
                        write_sequence( temp_seq , seq_filename , seqformat ,
                            os.path.split( root_name )[-1] + ' chain ' + chain ,
                             description , protein_alphabet , seqformat_extension_map )

                # write out the model+chain
                writer.set_structure( temp )
                print 'Writing a copy of model ' + str( model + 1 ) + ' chain(s) ' + complx + ' to ' + new_dir + '.pdb'
                writer.save( new_dir + '/' + dir_name + '.pdb' )#, selection )
                
                # also write a cleaned PDB file, onlt ATOM lines
                clean_ATOM_lines_from_pdb( new_dir + '/' + dir_name + '.pdb' )
                
                # also write any alternate conformations
                clean_alternate_conformations_from_pdb( new_dir + '/' + dir_name + '.pdb' )
            
                # also get specific HETATMs...this is getting bulky...
                clean_HETATM_lines_from_pdb( new_dir + '/' + dir_name + '.pdb' )
                
                # no need to clean DNA
    else:
        # only 1 model AND only 1 chain
        # still write it please :)
        model = 0
        chain = complexes[0]
        
        # may seem silly, but this edge case will prevent needless re-parsing

        # get the chain sequence
        seqid = os.path.split( root_name )[-1] + '_chain_' + complexes[0]

        extract_protein_sequence_from_pdb( struct , True ,
                seqid + '.' + seqformat_extension_map[seqformat] , seqid , '' ,
                model , chain , True ,
                seqformat = seqformat , alphabet = protein_alphabet , seqformat_extension_map = seqformat_extension_map )
    
    # debug summary...
    temp = os.listdir( os.getcwd() )
    temp.sort()
    print 'New Files in the ' + root_name + ' directory :\n' + '\n'.join( ['\t'+ i for i in temp] )
    
    # return back one directoy
    os.chdir( original_dir )    # yeah...its hacky
    
    return best_guess
    
################################################################################
# HEADER STUFF

# extract header text
def extract_header_from_pdb( pdb_filename , header_filename = 'header.txt' ):
    # write the header (?)
    # get the header
    f = open( pdb_filename , 'r' )
    header = ''
    while True:    # should error from f.next() if improper input...
        # next line
        line = f.next()
        
        # exit condition
        if 'ATOM' == line[:4] or 'MODEL' == line[:5] or 'HETATM' == line[:6]:
            break
        
        header += line
    f.close()
    
    # write the header
    if header_filename:
        print 'Writing a copy of the header lines to the file ' + header_filename
        f = open( header_filename , 'w' )
        f.write( header )
        f.close()

    return header

# return any predicted shain pairs
def extract_duplicate_chains_from_pdb_header( pdb_filename ):
    # load the raw data
    f = open( pdb_filename , 'r' )
    complexes = []
    keep_going = True
    while keep_going:
        # next line
        line = f.next()

        # ...think about this...
        # check if chain info, extract the matching subunits
        if line[:6] == 'COMPND' and 'CHAIN:' in line:
            duplicate = line.split( 'CHAIN: ' )[-1].replace( ';' , '' ).strip().split( ', ' )    # ignore ";\n"
            if len( duplicate ) > 1:
                complexes.append( duplicate )
        # stop condition
        elif not ('HEADER' in line or 'TITLE' in line or 'COMPND' in line or 'CAVEAT' in line):
            keep_going = False
    f.close()

    # convert complexes
    if complexes:
        if not sum( [len( c ) - len( complexes[0] ) for c in complexes] ):
            # all are the same length
            complexes = [''.join( [c[i] for c in complexes] ) for i in xrange( len( complexes[0] ) )]
        else:
            # uh oh...
            # could be all should be unique...which puts us in exception land anyway
            # assume that last listed are aberrantly unpaired
            lowest = min( [len( c ) for c in complexes] )
            temp = [''.join( [c[i] for c in complexes] ) for i in xrange( lowest )]
            for c in complexes:
                temp += c[lowest:]
            complexes = temp
    
    return complexes

# return number of models, scanned from header
def extract_number_of_models_from_pdb_header( pdb_filename ):
    # get the number of models
    f = open( pdb_filename , 'r' )
    models = 1
    keep_going = True
    while keep_going:
        # next line
        line = f.next()
        
        # check for models
        if line[:6] == 'NUMMDL':
            models = int( line.replace( 'NUMMDL' , '' ).strip() )
            keep_going = False
        elif line[:4] == 'ATOM':
            keep_going = False
    f.close()
    
    return models

# return resolution, scanned from header
# other information? R-value? R-free?
# other places to extract the quality...?
def extract_resolution_information_from_pdb_header( pdb_filename ):
    # load it
    f = open( pdb_filename , 'r' )
    
    # ewww....should be a "for" loop that breaks...
    keep_going = True
    experimental_data = 'X-RAY DIFFRACTION'
    resolution = None
    while keep_going:
        # next line
        line = f.next()
        
        # check for models
        if line[:6] == 'EXPDTA':
#            print 'found exp data'
            experimental_data = line[6:].strip()
        elif line[:10] == 'REMARK   2':
            # check for NMR
#            print 'found remark'
#            print line
            if 'ANGSTROMS' in line:
#                print 'found resolution'
                resolution = float( line[23:].strip().split( 'ANGSTROMS' )[0].strip() )
                keep_going = False
        elif line[:4] == 'ATOM':
            keep_going = False
    f.close()
    
    return resolution , experimental_data

# return number of models, scanned from header
def extract_HETNAM_from_pdb_header( pdb_filename ):
    # get the number of models
    f = open( pdb_filename , 'r' )
    hetname_map = {}
    keep_going = True
    while keep_going:
        # next line
        line = f.next()
        
        # check for models
        if line[:6] == 'HETNAM':
            hetname = line[6:].strip().split( ' ' )
            hetkey = hetname[0]
            hetname = ''.join( [i + ' ' for i in hetname[1:]] )[:-1]
            
            hetname_map[hetkey] = hetname
        elif line[:4] == 'ATOM':
            keep_going = False
    f.close()
    
    return hetname_map

################################################################################
# DIVIDE AND JOIN

# split or join PDB files

# simple wrapper
def morph_atomName2element( atomName ):
    """
    Returns the element in  <atomName>
    
    raw PDB atomNames are supposed to have the element as the first character
    """
    element = atomName[:2].strip()
    
    # remove number characters
    for i in '0123456789':
        element = element.replace( i , '' )
    
    return element

# make sure a filename, Structure, or Model returns the Model of interest
# not tested recently...
def load_pdb( pdb , model = 0 ):
    """
    Returns the  <model>  of  <pdb>  if its a Structure object (or a filename)
    """
    # sort the input
    if isinstance( pdb , str ):
        # filename
        print 'Input filename ' + pdb + ', loading the structure now'
        parser = PDBParser( PERMISSIVE = 1 )
        pdb = parser.get_structure( 'temp' , pdb )

        # default to first one if empty...
        if not model:
            model = pdb.child_dict.keys()[0]

        print 'extracting the first model (' + str( model ) + ')'
        pdb = pdb[model]    # get the first model

    # tried doing this a prettier way...
    # check for specific methods and data types for clues...
    elif isinstance( pdb.child_dict.keys()[0] , int ):
        # its a Biopython structure
        # default to first one if empty...
        if not model:
            model = pdb.child_dict.keys()[0]

        print 'Input Biopython Structure, extracting the first model (' + str( model ) + ')'
        pdb = pdb[model]    # get the first model

    elif 'child_dict' in dir( pdb ):
        # ...could be any number of things...including what we want!
        # hooray! everything is okay
        None
    else:
        # not supported!
        raise IOError( 'That data structure is not currently supported...' )
    
    return pdb

# check the PDB for models and split into separate PDBs
def split_pdb_into_models( pdb_filename ):
    """
    Writes a single PDB file for every model in  <pdb_filename>
    uses the Biopython PDBParser and PDBIO
    """
    # make tools
    parser = PDBParser( PERMISSIVE = 1 )
    writer = PDBIO()

    pdb_filename = os.path.abspath( pdb_filename )
    root_name = get_root_filename( pdb_filename )
    struct = parser.get_structure( root_name , pdb_filename )

    # over the models
    for i in struct.child_dict.keys():
        # get just the model
        temp = Structure( 'temp' )
        temp.add( struct[i] )

        # write it
        writer.set_structure( temp )
        out_filename = root_name + '_model_' + str( i + 1 ) + '.pdb'
        print 'Model ' + str( i + 1 ) + ' written to ' + out_filename
        writer.save( out_filename )

# check the PDB for chains and split into separate PDBs
def split_pdb_into_chains( pdb_filename , model = 0 , export = False ):
    """
    Writes a single PDB file for every chain in  <pdb_filename>
    uses the Biopython PDBParser and PDBIO
    """
    # make tools
    parser = PDBParser( PERMISSIVE = 1 )
    writer = PDBIO()

    pdb_filename = os.path.abspath( pdb_filename )
    root_name = get_root_filename( pdb_filename )
    struct = parser.get_structure( root_name , pdb_filename )
    
    # assume there is only 1 model
    # over the chains
    chains = {}
    for i in struct[model].child_dict.keys():
        # get just the model
        temp = Structure( 'temp' )
        temp_mod = Model( 0 )
        temp_mod.add( struct[0][i] )
        temp.add( temp_mod )

        # write it
        writer.set_structure( temp )
        out_filename = root_name + '_chain_' + i + '.pdb'
#        chains.append( 'Chain ' + i + ' written to ' + out_filename )
        chains[i] = out_filename
        writer.save( out_filename )

    # debug output
    for i in chains.keys():
        print 'Chain ' + i + ' written to ' + chains[i]

    # optionally export
    if export:
        return chains

# add all files together in the provided order
# not tested recently...
def join_pdb_files( files , out_filename = '' ):
    """
    Combines the contents of all  <files>  and writes it out to  <out_filename>
    
    a very simple method
    """
    # default filename
    out_filename_provided = True
    if not out_filename:
        out_filename_provided = False
    
    text = ''
    for i in files:
        # open it
        f = open( i , 'r' )
        
        # add the text
        text += f.read()
        f.close()
    
        # check if the name should be added
        if not out_filename_provided:
            if '.' in i:
                out_filename += i[:i.find( '.' )]
            else:
                out_filename += i

    # write the bastard love child
    f = open( out_filename , 'w' )
    f.write( text )
    f.close()

# extract the chains from the PDB
# only considers ATOM lines, mainly for use with clean_nucleic_acid_lines_from_pdb
def extract_chains_from_pdb( pdb_filename , only = ['ATOM'] ):
    """
    Returns the chains found in  <pdb_filename>
    
    Only consider lines starting with  <only>
    """
    pdb_filename = os.path.abspath( pdb_filename )
    if os.path.exists( pdb_filename ):
        # load the data
        f = open( pdb_filename , 'r' )
        data = [i for i in f.xreadlines() if i[:6].strip() in only]
        f.close()
    
        # find unique chains
        chains = []
        for i in data:
            if not i[21] in chains:
                chains.append( i[21] )

        return chains
    else:
        print 'No such file or directory named ' + os.path.relpath( pdb_filename )
        return False
    
# extract the name mapping
def extract_numbering_map_from_pdb( pdb_filename , out_filename = '' , only = ['ATOM'] ):
    """
    Returns a map (dict) from residues in  <pdb_filename>  that are 1-indexed
    and a reverse map (dict)
    
    Only consider lines starting with  <only>
    
    Optionally write the results to  <out_filename>
    """
    pdb_filename = os.path.abspath( pdb_filename )
    if os.path.exists( pdb_filename ):
        # load the raw data
        f = open( pdb_filename , 'r' )
        d = [i for i in f.xreadlines() if i[:6].strip() in only]
        f.close()
    
        # extract dict of pairs
        pdb_map = {}
        reverse_map = {}
        count = 0
        text = ''
        for i in d:
            # basic info
            chain = i[21]
            resseq = i[22:26].strip()
            icode = i[26]    # the icode
            key = chain + resseq + icode
        
            if not key in pdb_map.keys():
                count += 1
                pdb_map[key] = count
                reverse_map[count] = key
            
                text += key + '\t' + str( count ) + '\n'

        # optionally write to file
        # no defaulting!
        if out_filename:
            # default filename
    #        f = open( get_root_filename( pdb_filename ) + '_PDB_numbering.txt' , 'w' )
    #        f.write( ''.join( [i +'\t'+ str( pdb_map[i] ) +'\n' for i in pdb_map.keys()] )  )
            print 'Writing the PDB numbering of ' + pdb_filename + ' to ' + out_filename
            f = open( out_filename , 'w' )
            f.write( text )
            f.close()
    
        return pdb_map , reverse_map
    else:
        print 'No such file or directory named ' + os.path.relpath( pdb_filename )
        return False

# extract a protein sequence from a PDB
# make this better? specify the chain?
# for now, only works if single chain...
def extract_protein_sequence_from_pdb( pdb , include_breaks = True ,
        out_filename = '' , seqid = '' , description = '' ,
        model = 0 , chain = 'A' , export = True ,
        seqformat = SEQFORMAT , alphabet = ProteinAlphabet ,
        seqformat_extension_map = SEQFORMAT_EXTENSION_MAP ):
    """
    Returns the protein sequences found in  <pdb>  in  <model>
    
    Optionally  <export>  the sequence
    Optionally write to  <out_filename>  with  <seqid>
    
    note: does NOT scan for protein chains, it only dumps out the full
        protein sequence in the PDB file
        individual chains can be extracted using process_pdb
    """
    # ensure pdb is proper, must be a model
    pdb = load_pdb( pdb , model )    # necessary model?

    # format the chain input
    if not isinstance( chain , list ):
        chain = [chain]
   
    # over the desired chains
    # ugh...this should all be rewritten...
    sequences = []
    for c in chain:
        # get it
        if include_breaks:
            # extract the sequence as a Biopython Seq object
            # convert the model into a Structure, for getting the sequence
            for_seq = Structure( 'temp' )
            # ...oh yeah...must be model 0 and up >:[
            temp_model = Model( 0 )    # hardcoded...
            for_seq.add( temp_model )
#            for ch in pdb.child_dict.keys():
            # copy it all over directly
            for_seq[0].add( pdb[c] )

            # gap regions makred as "|"
            seq_builder = PPBuilder()
            pp = seq_builder.build_peptides( for_seq )
            seq = Seq( '|'.join( [str( frag.get_sequence() ) for frag in pp] ) , alphabet )
#            for frag in pp:
#                seq += frag.get_sequence() + '|'    # already a Biopython Seq
            seqr = SeqRecord( seq )
            seqr.description = description + ' missing residues (gap regions as \"|\")'*( '|' in seq )
        else:
            # just iterate and extract!
            seq = Seq( ''.join( [three2one[i.resname] for i in pdb.get_residues() if i.resname in three2one.keys() and i.get_parent().id == c] ) , alphabet )
            seqr = SeqRecord( seq )
            seqr.description = description
    
        # prepare to write
    #    seq = seq[:-1]
    #    seqr.description = 'missing residues (gap regions as \"|\")'*( '|' in seq )    # no need if no gaps
        seqr.id = seqid
        sequences.append( seqr )

    # optionally write the sequence
    if out_filename:
        write_sequence( sequences , out_filename , seqformat , alphabet , seqformat_extension_map )
    
    # optionally export the sequence
    if export:
        return get_sequence( sequences )
#        return str( seq )

# extract and write a file from the PDB
def extract_nucleic_acid_sequences_from_pdb( pdb_filename , out_filename = '' , NA = NUCLEIC_SEQUENCE_LETTERS_MAP , DNA = NA_CODES['DNA'] , seqformat = SEQFORMAT , alphabet = DNAAlphabet , seqformat_extension_map = SEQFORMAT_EXTENSION_MAP ):
    """
    Returns the protein sequences found in  <pdb_filename>
    
    Only consider resNames in  <NA>
    
    Optionally write to  <out_filename>
    """
    pdb_filename = os.path.abspath( pdb_filename )
    if os.path.exists( pdb_filename ):
        # load the data
        f = open( pdb_filename , 'r' )
        d = f.readlines()
        f.close()
    
        # print about fails/assumptions
        print 'Extracting nucleic sequences from ' + os.path.relpath( pdb_filename ) + '\nFor visibility, this method assumes A LOT!\n1. nucleotides are identified by a unique resSeq codes (with a proper resName)\n2. sequences are identified by unique chain IDs\n3. RNA is the default state\n4. DNA is identified by \"DG\" (etc.) OR \"T\" resi codes\n4. All sequences are continuous\n6. All sequences are recorded 5\' -> 3\' (and written to file in this order)'
    
        # check for nucleic lines - No, do while parsing
        # extract sequence
        NA_keys = NA.keys()
        DNA_keys = DNA.keys()
#        molecule = 'RNA'
        molecule_types = {}
        sequences = {}
        last = None
        for line in d:
            # must have C1 and a nucleic resi code to be considered a nucleotide
            resname = line[17:20].strip()
            resseq = line[22:27].strip()    # resseq
            if (line[:5] == 'ATOM ' or line[:4] == 'TER ') and resname in NA_keys:# and line[13:16].strip() == 'C1\'':
                # only novel lines
                if resseq == last:
                    continue
                last = resseq    # if the remainder will execute...
            
                # check for DNA
                chain = line[21]
                if [True for i in DNA_keys if i in resname]:
                    # its DNA
                    molecule_types[chain] = 'DNA'
                    # consider the whole chain DNA if ANY of the exclusive codes are present
                    # sometimes DNA is abbreviated without the "d" to designate "deoxy"
            
                # remember the letter
                if chain in sequences.keys():
                    # add the letter
                    sequences[chain] += NA[resname]    # map the code
                else:
                    # create it as well
                    sequences[chain] = NA[resname]
                    molecule_types[chain] = 'RNA'    # default

        # default out name
        root_filename = get_root_filename( pdb_filename )
        if not out_filename:
            out_filename = root_filename

        # write the sequences
        for chain in sequences.keys():
            # verify its not just a nucleotide
            seq = sequences[chain]
            if len( seq ) > 1:
                # determine the molecule type
                # record a proprt id
                seqr = SeqRecord( Seq( seq , alphabet ) )    # even if RNA (?)
                seqr.id = os.path.split( root_filename )[-1] + '_chain_' + chain
                seqr.description = molecule_types[chain]
            
                # oh yeah, write it, prints out by itself
                out_filename = seqr.id + '.' + seqformat_extension_map[seqformat]
                write_sequence( seqr , out_filename , seqformat , alphabet , seqformat_extension_map )

        return sequences , molecule_types    # empty dict will evaluate as false
    else:
        print 'No such file or directory named ' + os.path.relpath( pdb_filename )
        return False

################################################################################
# CLEANING METHODS

# HERE !!!

# a dirty input produces a cleaned output file :)
# default behavior is to produce output

# removes non ATOM lines from  <pdb_file>  and writes to  <out_file>
def clean_ATOM_lines_from_pdb( pdb_filename , out_filename = '' , HETATM_include = [] , excluded_atoms = ['CN'] , accepted_fields = ['ATOM ' , 'TER '] ):
    """
    Writes all lines in the PDB file  <pdb_filename>  beginning with "ATOM" or
    "TER" into  <out_filename>  (defaults to  <pdb_file>.clean.pdb)

    Optionally include HETATM lines with resNames in  <HETATM_include>

    Returns True if successful
    
    ...pretty much the same as:
    
        grep "ATOM" pdb_filename > out_filename

    example:
        clean_non_ATOM('1YY9.pdb')
    See also:
        Pose
        Pose.dump_pdb
        pose_from_pdb
        pose_from_rcsb
    """
    # get the file rootname
    pdb_filename = os.path.abspath( pdb_filename )
    root_filename = get_root_filename( pdb_filename )
    if not root_filename:    # in case it is improper, no "."
        root_filename = pdb_filename
    
    # an optional argument for PDB files not ending in .pdb
#    if not edit:
#        edit = 255

    # if the file exists
    if os.path.exists( pdb_filename ):
        # find all ATOM and TER lines
        f = open( pdb_filename , 'r' )
        data = f.readlines()
        f.close()
        good = []
        for i in data:
            if [True for j in accepted_fields if i[:len( j )] == j]:
#            if i[:5] == 'ATOM ' or i[:4] == 'TER ':
                # add your preference rules for ligands, DNA, water, etc.
                # check for excluded atoms
                if i[12:16].strip() in excluded_atoms:
                    # skip it, do not add to the list
                    continue
                
                good.append( i )
            elif i[:6] == 'HETATM' and i[17:20] in HETATM_include:
                # save for later, more processsing
                good.append( i )

        # stop condition
        if not good:
            # tell the user and exit
            print 'No ATOM or HETATM lines in ' + os.path.relpath( pdb_filename )
            return False

        # default output file to  <pdb_filename>.clean.pdb
        if not out_filename:
            out_filename = root_filename + '.clean.pdb'

        # write the found lines
        print 'if the file ' + os.path.relpath( out_filename ) + ' already exists, it will be overwritten!'
        f = open( out_filename , 'w' )
        f.writelines( good )
        f.close()

        print 'PDB file ' + os.path.relpath( pdb_filename ) + ' successfully cleaned, non-ATOM lines removed\nclean data written to ' + os.path.relpath( out_filename )
        return out_filename

    else:
        print 'No such file or directory named ' + os.path.relpath( pdb_filename )
        return False

# if you would prefer a simpler call using grep, it looks something like this
#    os.system("grep \"ATOM\" %s.pdb > %s.clean.pdb"%(pdb_file[:edit],pdb_file[:edit]))

# split the ATOM lines, only look for DNA lines
def clean_nucleic_acid_lines_from_pdb( pdb_filename , out_filename = '' , NA = NUCLEIC_SEQUENCE_LETTERS_MAP.keys() ):
    """
    Scan  <pdb_filename>  for any nucleic acid lines and writes these to
    <out_filename>
    
    defines nucleic acid resNames (three letter codes) as those with
    stripped codes in  <NA>
    
    default definition of nucleic acid resNames can be adjusted in settings.py
    """
    # get the file rootname
    pdb_filename = os.path.abspath( pdb_filename )
    root_filename = get_root_filename( pdb_filename )
    if not root_filename:    # in case it is improper, no "."
        root_filename = pdb_filename
    
    # if the file exists
    if os.path.exists( pdb_filename ):
        # find all ATOM and TER lines
        f = open( pdb_filename , 'r' )
        data = f.readlines()
        f.close()

        good = []
        for i in data:
            if (i[:5] == 'ATOM ' or i[:4] == 'TER ') and i[17:20].strip() in NA:
                # add your preference rules for ligands, DNA, water, etc.
                good.append( i )

        # stop condition
        if not good:
            # tell the user and exit
            print 'No nucleic acid lines in ' + os.path.relpath( pdb_filename )
            return False

        # default output file to  <pdb_filename>.clean.pdb
        if not out_filename:
            out_filename = root_filename + '.nucleic.pdb'

        # write the found lines
        print 'if the file ' + os.path.relpath( out_filename ) + ' already exists, it will be overwritten!'
        f = open( out_filename , 'w' )
        f.writelines( good )
        f.close()
        
        print 'PDB file ' + os.path.relpath( pdb_filename ) + ' successfully cleaned, DNA/RNA lines extracted\nclean data written to ' + os.path.relpath( out_filename )
        return out_filename

    else:
        print 'No such file or directory named '+ os.path.relpath( pdb_filename )
        return False

# split the ATOM lines, only look for not RNA/DNA lines
def clean_protein_lines_from_pdb( pdb_filename , out_filename = '' , NA = NUCLEIC_SEQUENCE_LETTERS_MAP.keys() ):
    """
    Scan  <pdb_filename>  for any nucleic acid lines and writes all "ATOM" lines
    that are NOt nucleic acids to  <out_filename>
    
    defines nucleic acid resNames (three letter codes) as those with
    stripped codes in  <NA>
    
    default definition of nucleic acid resNames can be adjusted in settings.py
    """
    # get the file rootname
    pdb_filename = os.path.abspath( pdb_filename )
    root_filename = get_root_filename( pdb_filename )
    if not root_filename:    # in case it is improper, no "."
        root_filename = pdb_filename
    
    # if the file exists
    if os.path.exists( pdb_filename ):
        # find all ATOM and TER lines
        f = open( pdb_filename , 'r' )
        data = f.readlines()
        f.close()

        good = []
        for i in data:
            if (i[:5] == 'ATOM ' or i[:4] == 'TER ') and not i[17:20].strip() in NA:
                # add your preference rules for ligands, DNA, water, etc.
                good.append( i )

        # stop condition
        if not good:
            # tell the user and exit
            print 'No protein lines in ' + os.path.relpath( pdb_filename )
            return False

        # default output file to  <pdb_filename>.clean.pdb
        if not out_filename:
            out_filename = root_filename + '.protein.pdb'

        # write the found lines
        print 'if the file ' + os.path.relpath( out_filename ) + ' already exists, it will be overwritten!'
        f = open( out_filename , 'w' )
        f.writelines( good )
        f.close()
        print 'PDB file ' + os.path.relpath( pdb_filename ) + ' successfully cleaned, protein lines extracted\nclean data written to ' + os.path.relpath( out_filename )
        return True

    else:
        print 'No such file or directory named '+ os.path.relpath( pdb_filename )
        return False

# scan for HETATMs, rewrite without all these lines, record specific ones
def clean_HETATM_lines_from_pdb( pdb_filename , out_filename = '' , only = '' , write_unique = True ):
    """
    Writes all lines in the PDB file  <pdb_filename>  beginning with "HETATM"
    into  <out_filename>  (defaults to  <pdb_filename>.hetatm.pdb)
    Optionally write PDB files for all unique residue type codes in the HETATM
    lines if  <write_unique>  is True (default True)
    
    OR
    
    Writes all lines in the PDB file  <pdb_filename>  beginning with "HETATM"
    AND with the resName  <only>

    Returns True if successful
    """
    # get the file rootname
    pdb_filename = os.path.abspath( pdb_filename )
    root_filename = get_root_filename( pdb_filename )
    if not root_filename:    # in case it is improper, no "."
        root_filename = pdb_filename

    # if the file exists
    if os.path.exists( pdb_filename ):
        # find all HETATM
        f = open( pdb_filename , 'r' )
        data = f.readlines()
        f.close()
        good = []
        unique = []
        for i in data:
            resn = i[17:20].strip()
            if i[:6] == 'HETATM' and (not only or resn in only):
                # save for later, more processsing
                good.append( i )
                
                # look for unique resn names
                if not only and not resn in unique:
                    unique.append( resn )
        
        # stop condition
        if not good:
            # tell the user and exit
            print 'No HETATM lines in ' + os.path.relpath( pdb_filename )
            return False
        
        # default output file to  <pdb_filename>.clean.pdb
        if not out_filename:        
            if not only:
                out_filename = root_filename + '.hetatm.pdb'
            elif only in WATER_CONVERSION.keys():    # just waters...
                out_filename = root_filename.replace( '.hetatm' , '' ) + '.waters.pdb'
            else:
                # its anything else, name based on the code
                out_filename = root_filename.replace( '.hetatm' , '' ) + '.' + only + '.pdb'

        # write the found lines
        print 'if the file ' + os.path.relpath( out_filename ) + ' already exists, it will be overwritten!'
        f = open( out_filename , 'w' )
        f.writelines( good )
        f.close()
        
        # change this!
        if not only:
            print 'PDB ' + os.path.relpath( pdb_filename ) + ' successfully cleaned, non-HETATM lines removed\nclean data written to ' + os.path.relpath( out_filename )
        else:
            print 'All ' + only + ' lines in PDB file ' + os.path.relpath( pdb_filename ) + ' written to ' + os.path.relpath( out_filename )

        # optionally redo for all unique members
        if not only and write_unique:
            if len( unique ) > 1:
                # do them all
#                for resn in unique:
#                    clean_HETATM_lines_from_pdb( out_filename , '' , resn )
                unique_filenames = [clean_HETATM_lines_from_pdb( out_filename , '' , resn ) for resn in unique]
                return out_filename , unique_filenames
            else:
                # only 1 HETATM type...
                unique = unique[0]
                print 'Only 1 type of HETATM found, ' + unique

                if unique in WATER_CONVERSION.keys():
                    unique = 'waters'

#                print 'Renaming ' + root_filename + '.hetatm.pdb to ' + root_filename + '.' + unique + '.pdb'
#                shutil.move( root_filename + '.hetatm.pdb' , root_filename + '.' + unique + '.pdb' )
                temp = root_filename + '.' + unique + '.pdb'
                print 'Renaming ' + os.path.relpath( out_filename ) + ' to ' + os.path.relpath( temp )
                shutil.move( out_filename , temp )
                out_filename = temp

        return out_filename
    else:
        print 'No such file or directory named ' + os.path.relpath( pdb_filename )
        return False

# scan for alternate location fields
def clean_alternate_conformations_from_pdb( pdb_filename , remove_identifier = True ):
    """
    Writes PDB files for each of the alternate conformations found in
    <pdb_filename>
    """
    # get the file rootname
    pdb_filename = os.path.abspath( pdb_filename )
    root_filename = get_root_filename( pdb_filename )
    if not root_filename:    # in case it is improper, no "."
        root_filename = pdb_filename
    
    # verify it exists
    if not os.path.exists( pdb_filename ):
        # for pipelines etc.
        print 'No such file or directory named ' + os.path.relpath( pdb_filename )
        return False
    
    # find all alternate conformations
    f = open( pdb_filename , 'r' )
    lines = f.readlines()
    f.close()
    
    # for storage
    non_alternating = ['']
    alternate_conformations = []
    last_line_alternate = False
    index = 0
    alternate_index = -1
    conformation_names = []
    resis = set()
    for i in lines:
        # skip non ATOM lines...fix this later to support header?
        if not i[:6].strip() in ['ATOM' , 'HETATM']:
            last_line_alternate = False
            continue
        
        # sort it
        if i[16].strip():
            conformation = i[16]
            resis.add( i[21] +':'+ i[22:27].strip() )
            
            # optionally remove the alternate conformation identifier
            if remove_identifier:
                i = i[:16] + ' ' + i[17:]
            
            # did we just transition into an alt conf region?
            if last_line_alternate:
                # still in the same region
                if not conformation in alternate_conformations[alternate_index].keys():
                    alternate_conformations[alternate_index][conformation] = i
                    if not conformation in conformation_names:
                        conformation_names.append( conformation )
                else:
                    alternate_conformations[alternate_index][conformation] += i
            else:
                # in a new region
#                if alternate_conformations:
#                    conformation_names = list( set( conformation_names + alternations_conformations[-1].keys() ) )
#                    number_of_conformations = max( number_of_conformations , len( alternate_conformations[-1].keys() ) )
                alternate_index += 1
                alternate_conformations.append( {conformation : i} )
                if not conformation in conformation_names:
                    conformation_names.append( conformation )
            
            last_line_alternate = True
        else:
            # did we just transition into an alt conf region?
            if last_line_alternate:
                # entered a new region
                index += 1
                non_alternating.append( i )
            else:
                # in the same region
                non_alternating[index] += i
                
            last_line_alternate = False
        
    # exit condition
    conformation_names.sort()    # intuitive order...
    if not conformation_names:
        print 'No alternate conformations detected (17th column)'
        return False
    else:
        print 'found ' + str( len( conformation_names ) ) + ' alternate conformations: ' + ', '.join( conformation_names )
        print 'alternate locations found for residues: ' + ', '.join( list( resis ) )
    
#    print index , alternate_index , number_of_conformations
    
    # write out the alternate conformations
    conformation_filenames = []
    for i in conformation_names:
        # make a text by building from fragments
        text = ''
        for j in xrange( len( non_alternating ) - 2 ):
            text += non_alternating[j]
            if i in alternate_conformations[j].keys():
                text += alternate_conformations[j][i]
            else:
                # default to the "first" alt conf ID
                key = 0
                while not conformation_names[key] in alternate_conformations[j].keys():
                    key += 1
                key = conformation_names[key]
                text += alternate_conformations[j][key]

        # add edge case
        text += non_alternating[-1]
        
        # write the file
        out_filename = root_filename + '_conformation_' + i +'.pdb'
        print 'writing conformation ' + i + ' out to ' + os.path.relpath( out_filename ) + ' ...'
        f = open( out_filename , 'w' )
        f.write( text )
        f.close()
        
        conformation_filenames.append( out_filename )

    return conformation_filenames



    

################################################################################
# CONVERTERS

# rewrite the hetatm lines in the pdb file
def convert_pdb_resnames_to_ATOM_lines( hetatm_pdb_filename , out_filename = '' , report_filename = '' , conversion = three2three ):
    """
    Rewrites all HETATM lines in  <hetatm_pdb_filename>  found as keys in
    the dict  <conversion>  and replaces them with their values
    also rewrites the "HETATM" record as "ATOM  "
    
    used to convert HETATM lines that are proxies for amino acids    
    """
    hetatm_pdb_filename = os.path.abspath( hetatm_pdb_filename )
    
    # handle defaults
    if not out_filename:
        # override
        print 'no output filename provided, overwriting ' + hetatm_pdb_filename
        out_filename = hetatm_pdb_filename
    
    # make sure it exists
    if os.path.isfile( hetatm_pdb_filename ):
        # load in the lines
        f = open( hetatm_pdb_filename , 'r' )
        d = f.readlines()
        f.close()
    
        # change to the desired format
        converted = []
        for line in xrange( len( d ) ):
            record = d[line][:6].strip()
            resname = d[line][17:20].strip()
            # go ahead and just rewrite
            if record in ['ATOM' , 'HETATM'] and not resname in one2three.values() and resname in conversion.keys():
                new = conversion[resname]
                d[line] = d[line][:17] + new.rjust(3) + d[line][20:]

                # for records...
                temp = resname + ' lines converted to ' + new
                if not temp in converted:
                    converted.append( temp )
                
                # check the record...all to ATOM
                if record == 'HETATM':
                    d[line] = 'ATOM  '+ d[line][6:]
        
        # debug output
        if converted:
            converted = '\n'.join( converted )
            print converted
            if report_filename:
                print 'summary of converted lines written to ' + report_filename
                f = open( report_filename , 'w' )
                f.write( converted )
                f.close()
        
        # write it back
        f = open( out_filename , 'w' )
        f.writelines( d )
        f.close()
    else:
        print 'No such file named ' + os.path.relpath( hetatm_pdb_filename )
        return False


# useful?

# rewrite the water lines in the pdb file to the standard...from settings?
def convert_water_containing_pdb( hetatm_pdb_filename , conversion = WATER_CONVERSION ):
    """
    Rewrites all HETATM "water" lines in  <hetatm_pdb_filename>  to resNames
    based on  <conversion>
    
    adjust the definition of water (<look_for>) and what to switch to in
    settings.py
    
    not currently used...
    """
    hetatm_pdb_filename = os.path.abspath( hetatm_pdb_filename )
    if os.path.isfile( hetatm_pdb_filename ):
        # load in the lines
        f = open( hetatm_pdb_filename , 'r' )
        d = f.readlines()
        f.close()
    
        # change to the desired format
        for line in xrange( len( d ) ):
            resname = d[line][17:20]
            if resname.strip() in WATER_CONVERSION.keys():
                d[line] = d[line][:17] + WATER_CONVERSION[resname].rjust(3) + d[line][20:]
    
        # write it back...bad!
        f = open( hetatm_pdb_filename , 'w' )
        f.writelines( d )
        f.close()
    else:
        print 'No such file named ' + os.path.relpath( hetatm_pdb_filename )
        return False


# removes lines from  <pdb_file>  and writes to  <out_file>  ending in new
def clean_ATOM_non_new_lines_from_pdb( pdb_filename , out_filename = '' ):
    """
    Write all lines in the PDB file  <pdb_filename>  as long as the last three
    characters on the line aren't "new"
    
    used to clean Hydrogens added using Reduce
    """
    # get the file rootname
    pdb_filename = os.path.abspath( pdb_filename )
    root_filename = get_root_filename( pdb_filename )
    if not root_filename:    # in case it is improper, no "."
        root_filename = pdb_filename
    
    # an optional argument for PDB files not ending in .pdb
#    if not edit:
#        edit = 255

    # if the file exists
    if os.path.exists( pdb_filename ):
        # find all ATOM and TER lines
        f = open( pdb_filename , 'r' )
        data = f.readlines()
        f.close()
        good = []
        for i in data:
            if (i[:5] == 'ATOM ' or i[:4] == 'TER ') and not i.strip()[-3:] == 'new' and i[17:20] in one2three.values():
                good.append( i )

        # stop condition
        if not good:
            # tell the user and exit
            print 'No ATOM non-new lines in ' + os.path.relpath( pdb_filename )
            return False

        # default output file to  <pdb_filename>.clean.pdb
        if not out_filename:
            out_filename = root_filename + '.non_new.pdb'

        # write the found lines
        print 'if the file ' + os.path.relpath( out_filename ) + ' already exists, it will be overwritten!'
        f = open( out_filename , 'w' )
        f.writelines( good )
        f.close()

        print 'PDB file ' + os.path.relpath( pdb_filename ) + ' successfully cleaned, non-ATOM lines lacking \"new\" removed\nclean data written to ' + os.path.relpath( out_filename )
        return out_filename

    else:
        print 'No such file or directory named ' + os.path.relpath( pdb_filename )
        return False

################################################################################
# MAIN

if __name__ == '__main__':
    # parser object for managing input options
    parser = optparse.OptionParser()
    # essential data
    parser.add_option( '-i' , dest = 'pdb_filename' ,
        default = '' ,
        help = 'the pdb filename to process' )

    (options,args) = parser.parse_args()

    process_pdb( options.pdb_filename )


################################################################################
################################################################################
# UNFINISHED!!!

# scan for repeated chains and delete them, rewrite it
def clean_redundancy_from_pdb( in_filename , out_filename = '' ):
    """
    Not currently supported
    """
    print 'This is not currently supported sorry...\nIt should look for redundant copies...though it seems the best way to do this is to read directly from the header...but...even for the same sequence, the PDB file may have slightly different coordinates..so how to choose?\nUse process_pdb instead, a separate method is not supported because of this choice problem'

# rewrite a dna or rna pdb to be rosetta friendly
def convert_nucleic_acids_for_rosetta( nucleic_pdb_filename ):
    """
    Not currently supported
    """
    print '...still researching...for whatever reason, most DNA PDB coordinates are accepted in Rosetta (and thus do not crash PyRosetta) however, I cannot get RNA loading to work no matter what (!!??!)\nthey can be ~loaded by modifying the database, but this does not seem to actually do anything, although...make_pose_from_sequence can then make RNA polymers, generating a nonstandard ResidueSet also does not work...perhaps cause the lines are ATOM?'


"""
deprecated stuff...just in case...

#    f = open( pdb_filename , 'r' )
#    complexes = []
#    keep_going = True
#    while keep_going:
        # next line
#        line = f.next()

        # ...think about this...
        # check if chain info, extract the matching subunits
#        if 'CHAIN:' in line:
#            dupl = line.split( 'CHAIN: ' )[-1].replace( ';' , '' ).strip().split( ', ' )    # ignore ";\n"
#            if len( dupl ) > 1:
#                complexes.append( dupl )
        # stop condition
#        elif not ('HEADER' in line or 'TITLE' in line or 'COMPND' in line):
#            keep_going = False        
#    f.close()

    
    # convert complexes
#    if complexes:
#        if not sum( [len( c ) - len( complexes[0] ) for c in complexes] ):
            # all are the same length
#            complexes = [''.join( [c[i] for c in complexes] ) for i in xrange( len( complexes[0] ) )]
#        else:
            # uh oh...
            # could be all should be unique...which puts us in exception land anyway
            # assume that last listed are aberrantly unpaired
#            lowest = min( [len( c ) for c in complexes] )
#            temp = [''.join( [c[i] for c in complexes] ) for i in xrange( lowest )]
#            for c in complexes:
#                temp += c[lowest:]
#            complexes = temp

#                        shutil.copy( seq_filename , seqid2 )
#                        extract_protein_sequence_from_pdb( temp , '' + seqid + '.fa' , seqid , model , chain )
                
                # so...this Biopython object just doesn't work...
                # make a new selection
    #            selection = Select()
    #            selection.accept_model( i )
    #            for l in c:
    #                selection.accept_chain( l )


    # return the filename of the "best" PDB made for Rosetta
    # also return PDB numbering map?
#    if has_nucleic:
#        pdb_filename = root_name + '/' + pdb_filename
#    else:
#        pdb_filename = root_name + '/' + root_name + '.clean.pdb'

    # extract numbering of the best
#    pdb_map , reverse_map = extract_numbering_map_from_pdb( pdb_filename , pdb_filename[:-4] + '_numbering_map.txt' )
    # uh oh, bug...dont wanna fix now
    # until this is proper, leave this to Rosetta...

#    return pdb_filename #, pdb_map



#    if not chain:
#        chain = pdb.child_dict.keys()[0]
    
    # copy the chain

    # conver the model into a Structure, for getting the sequence
#    for_seq = Structure( 'temp' )
    # ...oh yeah...must be model 0 and up >:[
#    temp_model = Model( 0 )
#    for_seq.add( temp_model )
#    for ch in pdb.child_dict.keys():
        # copy it all over directly
#        for_seq[0].add( pdb[ch] )

    # extract the sequence as a Biopython Seq object
    # gap regions makred as "|"
#    seq_builder = PPBuilder()
#    pp = seq_builder.build_peptides( for_seq )
#    seq = Seq( '' , ProteinAlphabet )
#    for frag in pp:
#        seq += frag.get_sequence() + '|'    # already a Biopython Seq

    # or...just do this...

# from making sequences for subdirectories...
#                        temp_seq = SeqRecord( Seq( temp_seq , protein_alphabet ) )
#                        temp_seq.id = os.path.split( root_name )[-1] + ' chain ' + chain
#                        temp_seq.description = '(from model ' + str( model + 1 ) + ')'

"""


