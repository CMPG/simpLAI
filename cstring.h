#ifndef _MY_Cstring_HPP_
#define _MY_Cstring_HPP_


#ifdef _GCC_
#include <iostream>  //for << operator
#include <cstring>   //for functions strlen strcmp strchr
#else
#include <iostream.h>  //for << operator     
#include <cstring.h>   //for functions strlen strcmp strchr
#endif
#include <fstream>   //for read_token
#include <cstdio>    //for sprint

using namespace std;

const int NPOS = -1;

void my_strrev(char* s); //reverse the char* array  (not ANSI C++)

const int maxStringLength = 210000;

class my_string {

private:
    char* _data;
    int _del_inc;
    int _nb_char;
    int _tot_space;
    char* _reserve;

protected:
    int tot_space() {
        return _tot_space;
    }
public:

    void initialize() {
        _data = new char[10];
        _data[0] = '\0';
        _del_inc = 100;
        _tot_space = 10;
        _nb_char = 0;
        _reserve = NULL;
    }

    my_string() {initialize(); }

    my_string(const char& c) {
        _data = new char[10];
        _data[0] = c;
        _data[1] = '\0';
        _del_inc = 100;
        _tot_space = 10;
        _nb_char = 1;
        _reserve = NULL;
    }

    my_string(int i) {
        _data = new char[i + 2];
        _data[0] = '\0';
        _del_inc = 100;
        _tot_space = i + 2;
        _nb_char = 0;
        _reserve = NULL;
    }

    my_string(const my_string& x) {
        _data = new char[10];
        _data[0] = '\0';
        _del_inc = 100;
        _tot_space = 10;
        _nb_char = 0;
        if (this!=&x) *this = x;
        _reserve = NULL;
    }

    my_string(char* x) {
        _nb_char = strlen(x);
        _del_inc = 100;
        _data = new char[_nb_char + 2];
        strcpy(_data, x);
        _tot_space = _nb_char + 2;
        _reserve = NULL;
    }

    my_string(const char* x) {
        _nb_char = strlen(x);
        _del_inc = 100;
        _data = new char[_nb_char + 2];
        if (_nb_char) strcpy(_data, x);
        else _data[0]='\0';
        _tot_space = _nb_char + 2;
        _reserve = NULL;
    }

    ~my_string() {
        if (_data) delete[] _data;
        if (_reserve) delete[] _reserve;
    }

    char* c_str() const {
        return _data;
    }

    int length() const {
        return _nb_char;
    }
    void Lower();
    void Upper();
    void read_token(istream& is);
    //void read_token(strstream& is);
    my_string & operator=(const my_string& x);
    my_string & operator=(char* x);
    my_string & operator=(const char& x);
    my_string & operator=(int x);
//    my_string & operator=(float x);
    my_string & operator=(double x);
    my_string & operator+=(const my_string& x);
    my_string & operator+=(char* x);
    my_string & operator+=(char x); // new
    my_string & operator+=(int & x);

    char operator[](int pos) const {
        return _data[pos];
    } //new

    char& operator[](int pos) {
        return _data[pos];
    } //new
    friend my_string operator+(const my_string& x, const my_string& y);
    friend my_string operator+(const my_string& x, char* y);
    friend ostream & operator<<(ostream& os, const my_string& x);
    friend istream & operator>>(istream& is, my_string& s);

    int operator==(const my_string& x) const {
        return !(strcmp(_data, x.c_str()));
    }

    int operator==(char* x) const {
        return !(strcmp(_data, x));
    }

    int operator!=(const my_string& x) const {
        return (strcmp(_data, x.c_str()));
    }

    int operator!=(char* x) const {
        return (strcmp(_data, x));
    }

    int operator<(const my_string& x)const {
        return ((strcmp(_data, x.c_str())) < 0);
    } //new

    int operator>(const my_string& x)const {
        return ((strcmp(_data, x.c_str())) > 0);
    } //new
    //added since 21_3_97:
    void to_lower(); //converts my_string to lower case
    void to_upper(); //converts my_string to upper case
    void read_to_delim(istream& is, char delim); //read to delimiter
    void read_to_delim(istream& is); //read to end of line !!!
    //Loro_02_03_10
    void read_to_separator(istream& is, bool ignoreInitialSeps=false); //read to white char or tab

    my_string extract_sub_str(char delim); //extracts until the next occurence of delim or end of string //new 11_11_97
    my_string& remove(int pos, int n);
    int find_first_of(const char& s, int pos) const;
    int find_first_of(const char& s) const;
    int read_line(istream& is);
    int is_null();
    int contains(const char* pat) const;
    int contains(char pat) const;
    int contains(const my_string& s) const;

    void assign(const my_string& s) {
        *this = s;
    }
    int rfind(const my_string& s); //return the position of pattern s if not found return -1
    void rm_path(); //stef_28_3_98 removes path info from a file name

    void remove_blanks();
    void remove_repeated_seperators();
    my_string& remove_extension(); //Loro_16_01_01
    my_string& remove_path(); //Loro_18_02_15
    void remove_trailing_blanks(); //Loro_20_03_14
    
    void remove_comments(); //Loro_20_03_14
    void read_line_and_remove_comments(istream & is) {
       this->read_line(is);
       this->remove_repeated_seperators();
       this->remove_comments();
       this->remove_trailing_blanks();               
    }

    my_string extract_after(char delim); //returns the part of the string after the char delim
    //if delim does not exist, it returns the whole string
    my_string extract_after(int index);
    my_string extract_before(char delim); //returns the part of the string before the char delim
    //if delim does not exist, it returns the whole string
    my_string extract_before_doubleSlash(); // returns the part of the string before the double slash "//"
    my_string extract_before_hashtag(); // returns the part of the string before hashtag "#"
    my_string extract_before(int index);
    my_string extract_between(int first, int last);

};

#include "arrays.h"
typedef MY_TArrayAsVector <my_string> stringArray;

#endif







