/* The header file WvSpectr.h to declare the WvSpectr class and its related 
 * functions, including extracting spectra or one element from all cells.
 *                                                           JGLi02Jul2026  */

#ifndef WvSpectr_H
#define WvSpectr_H

#include <vector>
#include <cstddef>

class WvSpectr {

public:

    WvSpectr(size_t NSea, size_t NFrq, size_t NDrc);

    // Element access
    float& operator()(size_t i, size_t j);
    const float& operator()(size_t i, size_t j) const;

    // Spectral access
    std::vector<float>& operator[](size_t i);
    const std::vector<float>& operator[](size_t i) const;

    // Component operations
    std::vector<float> Compnt(size_t j);
    void setCompnt(size_t j, std::vector<float>& comp);

    // Inline size accessors
    size_t NSeapt() const { return NSe; }
    size_t NFrqcy() const { return NFr; }
    size_t NDrctn() const { return NDr; }
    size_t MSpctr() const { return NFr*NDr; }

    // Friend function to initialize spectral field.
    friend int Spectini( WvSpectr& Spc );
    // Friend function to integrate spctrum for significant wave height.
    friend std::vector<float> SWHeight( WvSpectr& Spc ); 

private:

    const size_t NSe;
    const size_t NFr;
    const size_t NDr;

    std::vector<std::vector<float>> spctrs;

};   // End of class WvSpectr declaration.

#endif

