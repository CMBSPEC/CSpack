//====================================================================================================================
//
// simple parser written by Boris Bolliet and Jens Chluba (Feb 2018). This was motivated by a parser from Class.
//
// Purpose: read parameters with given identifier from parameter file. It is assumed that the entry is written in the
//          form 'id-string = entry'. Lines starting with '#' and empty lines are omitted. For entries of the form
//          'id-string = entry entry2' the read functions currently only account for the first entry.
//
//====================================================================================================================
// 23.04.2018: overcame problem with variable names appearing after '#' with was not at beginning of line

//====================================================================================================================
// Standards
//====================================================================================================================
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include "parser.h"

using namespace std;

//====================================================================================================================
int parser_read_file(string filename, struct file_content &pfc, bool show_lines)
{
    ifstream ifile(filename.c_str());
    
    if(!ifile) {
        cerr << " parser_read_file :: Error opening parameter file " << filename << " Exiting. " << endl;
        exit(1);
    }

    pfc.filename=filename;
    
    do {
        string line="";
        getline(ifile, line);
        // drop comment lines
        if(line[0]!='#' && line.length()>1)
        {
            pfc.lines.push_back(line);
            // erase things after possible second '#'
            size_t pos=line.find('#');
            if(pos==string::npos) continue;
            pfc.lines.back().erase(pos, line.length());
        }
    }
    while(!ifile.eof());
    
    ifile.close();
    
    if(show_lines) for(int l=0; l<(int)pfc.lines.size(); l++) cout << pfc.lines[l] << endl;
    
    return 0;
}

//====================================================================================================================
int parser_free(struct file_content &pfc)
{
    pfc.lines.clear();
    return 0;
}

//====================================================================================================================
template <class T>
int parser_read_T(const struct file_content &pfc, string var_id, T &val, bool &found, bool show_entry)
{
    found = 0;
    
    for(int l=0; l<(int)pfc.lines.size(); l++)
    {
        string str=pfc.lines[l];
        size_t pos=str.find(var_id);
        size_t pos_found=pos;
        
        // continue if at end of string
        if(pos>=str.length() || pos_found>0) continue; // fixed problem with finding partial strings [JC]
     
        // if string is found, continue search from positions on
        pos=str.find("=", pos)+1;
        // check that this was really the callid and not a 'value'
        if(pos<=pos_found) continue;
        
        for(; pos<str.length(); pos++) if(str[pos]!=' ' || str[pos]=='#') break;
        
        // now at position of entry and need to convert it (everything after entry is omitted)
        istringstream iss(str.substr(pos));
        iss >> val;
        found = 1;
    }
    
    if(show_entry==1 && found==1) cout << var_id << " = " << val << endl;
    
    return 0;
}

int parser_read(const struct file_content &pfc, string var_id, int &val, bool &found, bool show_entry)
{ return parser_read_T(pfc, var_id, val, found, show_entry); }

int parser_read(const struct file_content &pfc, string var_id, double &val, bool &found, bool show_entry)
{ return parser_read_T(pfc, var_id, val, found, show_entry); }

int parser_read(const struct file_content &pfc, string var_id, string &val, bool &found, bool show_entry)
{ return parser_read_T(pfc, var_id, val, found, show_entry); }

//====================================================================================================================
int parser_read_fullstring(const struct file_content &pfc, string var_id, string &val,
                           bool &found, bool show_entry)
{
    found = 0;
    
    for(int l=0; l<(int)pfc.lines.size(); l++)
    {
        string str=pfc.lines[l];
        size_t pos=str.find(var_id);
        
        // continue if at end of string
        if(pos>=str.length()) continue;
     
        // if string is found, continue search from positions on
        pos=str.find("=", pos)+1;
        for(; pos<str.length(); pos++) if(str[pos]!=' ') break;
        
        size_t pos_start=pos;
        for(; pos<str.length(); pos++) if(str[pos]=='#') break;
        
        val = str.substr(pos_start, pos);
        
        found = 1;
    }
    
    if(show_entry==1 && found==1) cout << var_id << " = " << val << endl;
    
    return 0;
}

//====================================================================================================================
int parser_read(const struct file_content &pfc, string var_id,
                vector<double> &dvals,
                bool &found, bool show_entry)
{
    found = 0;
    double val;
    dvals.clear();
    
    for(int l=0; l<(int)pfc.lines.size(); l++)
    {
        string str=pfc.lines[l];
        size_t pos=str.find(var_id);
        size_t pos_found=pos;
        
        // continue if at end of string
        if(pos>=str.length()) continue;
     
        // if string is found, continue search from positions on
        pos=str.find("=", pos)+1;
        // check that this was really the callid and not a 'value'
        if(pos<=pos_found) continue;
        
        for(; pos<str.length(); pos++) if(str[pos]!=' ' || str[pos]=='#') break;
        
        // now at position of entry and need to convert it (everything after entry is omitted)
        istringstream iss(str.substr(pos));
        do
        {
            iss >> val;
            if(dvals.size()==0) dvals.push_back(val);
            else if(val!=dvals.back()) dvals.push_back(val);
        }
        while(!iss.eof());
        
        found = 1;
    }
    
    if(show_entry==1 && found==1)
    {
        cout << var_id << " = ";
        for(int k=0; k<(int)dvals.size(); k++)
            cout << dvals[k] << " ";
        cout << endl;
    }
    
    return 0;
}

//====================================================================================================================
int parser_read(const struct file_content &pfc, string var_id,
                vector<string> &svals,
                bool &found, bool show_entry)
{
    found = 0;
    string val;
    svals.clear();
    
    for(int l=0; l<(int)pfc.lines.size(); l++)
    {
        string str=pfc.lines[l];
        size_t pos=str.find(var_id);
        size_t pos_found=pos;
        
        // continue if at end of string
        if(pos>=str.length()) continue;
     
        // if string is found, continue search from positions on
        pos=str.find("=", pos)+1;
        // check that this was really the callid and not a 'value'
        if(pos<=pos_found) continue;
        
        for(; pos<str.length(); pos++) if(str[pos]!=' ' || str[pos]=='#') break;
        
        // now at position of entry and need to convert it (everything after entry is omitted)
        istringstream iss(str.substr(pos));
        do
        {
            iss >> val;
            if(svals.size()==0) svals.push_back(val);
            else if(val!=svals.back()) svals.push_back(val);
        }
        while(!iss.eof());
        
        found = 1;
    }
    
    if(show_entry==1 && found==1)
    {
        cout << var_id << " = ";
        for(int k=0; k<(int)svals.size(); k++)
            cout << svals[k] << " ";
        cout << endl;
    }
    
    return 0;
}

//====================================================================================================================
template <class T>
int parser_dual_read_T(const struct file_content &pfc, string var_id_1, string var_id_2,
                       T &val, bool &found_1, bool &found_2, bool show_entry)
{
    // Parser two similar parameters and pick the one that is actually found, else set an error.
    // If neither are found, then just set the found variable as
    parser_read(pfc, var_id_1, val, found_1, show_entry);
    parser_read(pfc, var_id_2, val, found_2, show_entry);

    if (found_1 && found_2)
    {
        string errmsg="Error! You have given "+var_id_1+" and "+var_id_2+". Please only supply one of these!";
        throw_error("parser_dual_read", errmsg);
    }
    return 0;
}

int parser_dual_read(const struct file_content &pfc, string var_id_1, string var_id_2,
                     int &val, bool &found_1, bool &found_2, bool show_entry)
{
    return parser_dual_read_T(pfc, var_id_1, var_id_2, val, found_1, found_2, show_entry);
}

int parser_dual_read(const struct file_content &pfc, string var_id_1, string var_id_2,
                     double &val, bool &found_1, bool &found_2, bool show_entry)
{
    return parser_dual_read_T(pfc, var_id_1, var_id_2, val, found_1, found_2, show_entry);
}

int parser_dual_read(const struct file_content &pfc, string var_id_1, string var_id_2,
                     string &val, bool &found_1, bool &found_2, bool show_entry)
{
    return parser_dual_read_T(pfc, var_id_1, var_id_2, val, found_1, found_2, show_entry);
}

//====================================================================================================================
//====================================================================================================================
