//
// Created by eason on 5/25/21.
//

#ifndef DNA_ENCODING_TRANSFORMATION_H
#define DNA_ENCODING_TRANSFORMATION_H

// transformation -- swap
// g_swap_granularity = how many nts as a unit to do swap
string swap(string strand) {
    if (g_swap_granularity==0) return strand;

    for (int i = 0; i < strand.size()-2*g_swap_granularity; i+=2*g_swap_granularity) {
        int first = i;
        int second = i+g_swap_granularity;

        string first_half=strand.substr(first,g_swap_granularity);
        string second_half=strand.substr(second,g_swap_granularity);
        // start point, len, const string
        strand.replace(first,g_swap_granularity,second_half);
        strand.replace(second,g_swap_granularity,first_half);
    }
    return strand;
}


// TODO: implement 3 mappings
string mapping(string strand) {
    string result;
    for (int i = 0; i < strand.size(); i++) {
        char nt = strand[i];
        /*mapping
        A -> A
        T -> T
        G -> C
        C -> G
        */
        /*if (nt == 'A') nt='T';
        else if(nt == 'T') nt='A';
        if (nt == 'C') nt='G';
        else if(nt == 'G') nt='C';*/
        result+=nt;
    }
    return result;
}
#endif //DNA_ENCODING_TRANSFORMATION_H
