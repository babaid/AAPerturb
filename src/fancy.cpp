#include "fancy.h"

#include<iostream>
#include<string>


// https://github.com/GregoryConrad/pBar/tree/master

ProgressBar::ProgressBar(int neededprogress): neededProgress(neededprogress){
}

void ProgressBar::update() {
    currentProgress += 1;
    amountOfFiller = (int)((currentProgress / neededProgress)*(double)pBarLength);
}

void ProgressBar::print(const std::string& message) {
    currUpdateVal %= pBarUpdater.length();
    std::cout << "\r" //Bring cursor to start of line
              << firstPartOfpBar; //Print out first part of pBar
    for (int a = 0; a < amountOfFiller; a++) { //Print out current progress
        std::cout << pBarFiller;
    }
    if(currentProgress!=neededProgress)
    {
        std::cout << pBarUpdater[currUpdateVal];
    }
    for (int b = 0; b < pBarLength - amountOfFiller; b++) { //Print out spaces
        std::cout << " ";
    }
    std::cout << lastPartOfpBar //Print out last part of progress bar
              << " (" << (int)(100*(currentProgress/neededProgress)) << "% | " << message << ")"
              //This just prints out the percent
              << std::flush;
    currUpdateVal += 1;
}
