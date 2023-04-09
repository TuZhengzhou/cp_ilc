#include <iostream>
#include <unordered_map>
#include <string> 
#include "structs.hpp"
using namespace std;

string linearEncode(string message) {
    unordered_map<char, int> charCount;
    for (char c : message) {
        charCount[c]++;
    }

    string code;
    for (char c : message) {
        if (charCount[c] == 0) continue; // already encoded
        code += to_string(charCount[c]); // add count to code
        code += c; // add character to code
        charCount[c] = 0; // mark as encoded
    }

    return code;
}

int main() {
    string message = "abbcccddddeeeee";
    string code = linearEncode(message);
    cout << "Original message: " << message << endl;
    cout << "Encoded message: " << code << endl;


    return 0;
}
