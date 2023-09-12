//
// Created by babaid on 09.09.23.
//
#include "molecules.h"


bool Atom::operator==(const Atom &atom) const {
    bool ser = this->serial == atom.serial;
    bool name = this->name == atom.name;
    bool rname = this->resName == atom.resName;
    bool el = this->element == atom.element;
    return ser && name && rname && el;
}