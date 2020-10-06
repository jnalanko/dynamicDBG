/**

Copyright 2017 Alan Kuhnle, Victoria Crawford

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <iostream>
#include <cmath>

/**
 * Implements an array of bits by storing integers
 */
class BitArray
{

private:
    // pointer to all our ints, which hold the bits
    uint64_t *ints;

    // the number of bits in an int
    uint64_t int_size;

    // the number of ints that ints points to
    uint64_t num_ints;

public:
    // Create a bit array of n bits
    BitArray(size_t num_bits)
    {

        allocate(num_bits);
    }

    BitArray()
    {
        //default constructor
    }

    void allocate(size_t num_bits)
    {
        // size of an uint64_t is given in bytes
        this->int_size = 8 * sizeof(uint64_t);

        this->num_ints = ceil((num_bits) / ((double)this->int_size));

        // create an array of that number of ints
        this->ints = new uint64_t[num_ints];

        // Make sure all ints are 0
        this->clearInts();
    }

    void clearInts()
    {

        for (uint64_t i = 0; i < this->num_ints; i++)
        {
            this->ints[i] = 0;
        }
    }

    // Get which int a certain bit is in
    uint64_t int_num(size_t bit_num)
    {

        return bit_num / this->int_size;
    }

    // Get what index of its int a certain bit is in
    uint64_t int_index(size_t i)
    {

        return i % this->int_size;
    }

    // get the ith bit
    bool get(size_t i)
    {

        // What int this bit_num is in
        uint64_t bit_int = this->ints[this->int_num(i)];
        // What spot in that int this bit_num is
        uint64_t bit_ind = this->int_index(i);

        // get this bit in the last spot of an int
        uint64_t num = bit_int >> (this->int_size - bit_ind - 1);
        // Then get higher order bits off
        num &= 1;

        return (bool)num;
    }

    // set the ith bit to v
    void set(size_t i, bool v)
    {

        // which int in our int array
        uint64_t int_num = this->int_num(i);

        // what index from the right
        uint64_t int_ind_right = this->int_size - this->int_index(i) - 1;

        // Set that spot
        uint64_t op = 1;
        op = op << int_ind_right;

        if (v)
        {
            // set to 1
            this->ints[int_num] |= op;
        }
        else
        {
            // set to 0
            op = ~op;
            this->ints[int_num] &= op;
        }
    }

    // Get the amount of bits this is taking up
    size_t total_bit_size()
    {

        return static_cast<size_t>(this->num_ints) * this->int_size;
    }

    // Print a string of bits to screen
    void print()
    {

        uint64_t index = 0;

        for (int64_t i = 0; i < this->num_ints; i++)
        {

            for (int64_t j = 0; j < this->int_size; j++)
            {
                std::cout << this->get(index);
                index++;
            }

            std::cout << std::endl;
        }
    }

    /* writes BitArray to a binary file stream */
    void save(ostream &of)
    {
        of.write((char *)(&num_ints), sizeof(uint64_t));
        of.write((char *)(&int_size), sizeof(uint64_t));

        uint64_t *ptr = ints;
        for (uint64_t i = 0; i < num_ints; ++i)
        {
            of.write((char *)(ptr++), sizeof(uint64_t));
        }
    }

    void load(istream &of)
    {
        of.read((char *)(&num_ints), sizeof(uint64_t));
        of.read((char *)(&int_size), sizeof(uint64_t));
        //    cerr << "Bitarray int_size " << int_size << endl;
        //    cerr << "Bitarray num_ints " << num_ints << endl;
        ints = new uint64_t[num_ints];
        uint64_t *ptr = ints;
        for (uint64_t i = 0; i < num_ints; ++i)
        {
            of.read(((char *)(ptr++)),
                    sizeof(uint64_t));
        }
    }

    ~BitArray()
    {

        delete[] this->ints;
    }
};
