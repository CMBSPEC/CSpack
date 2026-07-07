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

#ifndef PARSER_H
#define PARSER_H

#include <string>
#include <vector>
#include "routines.h"

using namespace std;

struct file_content
{
    vector<string> lines;
    string filename;
};

//====================================================================================================================
int parser_read_file(string filename, struct file_content &pfc, bool show_lines=0);
int parser_free(struct file_content &pfc);

int parser_read(const struct file_content &pfc, string var_id, int &val, bool &found, bool show_entry=0);
int parser_read(const struct file_content &pfc, string var_id, double &val, bool &found, bool show_entry=0);
int parser_read(const struct file_content &pfc, string var_id, string &val, bool &found, bool show_entry=0);

int parser_read_fullstring(const struct file_content &pfc, string var_id, string &val,
                           bool &found, bool show_entry=0);

int parser_read(const struct file_content &pfc, string var_id, vector<double> &dvals,
                bool &found, bool show_entry=0);

int parser_read(const struct file_content &pfc, string var_id, vector<string> &svals,
                bool &found, bool show_entry=0);

int parser_dual_read(const struct file_content &pfc, string var_id_1, string var_id_2,
                     int &val, bool &found_1, bool &found_2, bool show_entry);
int parser_dual_read(const struct file_content &pfc, string var_id_1, string var_id_2,
                     double &val, bool &found_1, bool &found_2, bool show_entry);
int parser_dual_read(const struct file_content &pfc, string var_id_1, string var_id_2,
                     string &val, bool &found_1, bool &found_2, bool show_entry);


#endif

//====================================================================================================================
//====================================================================================================================
