//
// Created by babaid on 14.09.23.
//
// This header contains all the fancy stuff which is fun and helpful to use. Like a progressbar.
// Not really important, but nice to look at.
// Fancy formatted stuff.
// Most of it not from me, just edited it to suit my needs.

#ifndef AA_PERTURB_FANCY_H
#define AA_PERTURB_FANCY_H
#include<string>
class ProgressBar {
public:
    ProgressBar(int neededprogress);
    void update();
    void print(std::string);
    std::string firstPartOfpBar = "[",
                lastPartOfpBar = "]",
                pBarFiller = "#",
                pBarUpdater = "/-\\|";
private:
    int amountOfFiller,
        pBarLength = 50,
        currUpdateVal = 0;
    double  currentProgress = 0,
            neededProgress;
};

#endif