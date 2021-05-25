//
// Created by eason on 5/25/21.
//

#ifndef DNA_ENCODING_TOOL_H
#define DNA_ENCODING_TOOL_H
bool isDir(string dir)
{
    struct stat fileInfo;
    stat(dir.c_str(), &fileInfo);
    if (S_ISDIR(fileInfo.st_mode)) {
        return true;
    } else {
        return false;
    }
}


// Function to convert a decimal
// number to a ternary number
string convertToTernary(int N)
{
    // Base case
    if (N == 0)
        return "0";

    // Finding the remainder
    // when N is divided by 3
    int x = N % 3;
    N /= 3;

    // Recursive function to
    // call the function for
    // the integer division
    // of the value N/3
    string result = convertToTernary(N);

    return result+to_string(x);
}

// Function to convert the decimal to ternary
string convert(int Decimal)
{
    // If the number is greater
    // than 0, compute the
    // ternary representation
    // of the number
    if (Decimal != 0) {
        return convertToTernary(Decimal);
    }
    else
        return "000000";
}
#endif //DNA_ENCODING_TOOL_H
