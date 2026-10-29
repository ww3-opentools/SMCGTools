/* The WvSpectr.cpp file defines the full/sub domain wave spectral class and 
 * some access/reset functions for a full spectrum or a compenent of all sea
 * points within the domain.                                JGLi02Jul2026   */

#include "WvSpectr.h"
#include <stdexcept>

// Constructor

WvSpectr::WvSpectr(size_t NSea, size_t NFrq, size_t NDrc)
    : NSe(NSea), NFr(NFrq), NDr(NDrc),
      spctrs(NSea, std::vector<float>(NFrq*NDrc)) {}

// Element access operator ()
float& WvSpectr::operator()(size_t iSe, size_t jSp) {
    return spctrs.at(iSe).at(jSp);
}

const float& WvSpectr::operator()(size_t iSe, size_t jSp) const {
    return spctrs.at(iSe).at(jSp);
}

// Spectral access at sea point iSe.
std::vector<float>& WvSpectr::operator[](size_t iSe) {
    return spctrs.at(iSe);
}
const std::vector<float>& WvSpectr::operator[](size_t iSe) const {
    return spctrs.at(iSe);
}

// Copy component operations
std::vector<float> WvSpectr::Compnt(size_t jSp) {
    size_t MSp = NFr*NDr;
    if (jSp >= MSp) {
        throw std::out_of_range("Spectral index out of range");
    }

    std::vector<float> comp(NSe);

    for (size_t i = 0; i < NSe; ++i) {
        comp[i] = spctrs[i][jSp];
    }

    return comp;
}

// Reset component operations.
void WvSpectr::setCompnt(size_t jSp, std::vector<float>& comp) {
    size_t MSp = NFr*NDr;
    if (jSp >= MSp) {
        throw std::out_of_range("Component index out of range");
    }

    if (comp.size() < NSe) {
        throw std::invalid_argument("Input component size mismatch");
    }

    for (size_t i = 0; i < NSe; ++i) {
        spctrs[i][jSp] = comp[i];
    }
}

