#ifndef MESHIO_H
#define MESHIO_H

#include <string>
#include <iostream>
#include <fstream>
#include <istream>
#include "domain.h"

namespace Oceane {

class MeshIO
{
public:
    MeshIO()=default;
    void strip_ws(std::string & upper)
    {
        upper.erase(std::remove_if(upper.begin(), upper.end(), isspace), upper.end());
    }

    std::string parse_label(std::string line, std::string label_name)
    {
      strip_ws(line);

      // Do all string comparisons in upper-case
      std::string
        upper_line(line),
        upper_label_name(label_name);
      std::transform(upper_line.begin(), upper_line.end(), upper_line.begin(), ::toupper);
      std::transform(upper_label_name.begin(), upper_label_name.end(), upper_label_name.begin(), ::toupper);

      // Get index of start of "label="
      size_t label_index = upper_line.find(upper_label_name + "=");

      if (label_index != std::string::npos)
        {
          // Location of the first comma following "label="
          size_t comma_index = upper_line.find(",", label_index);

          // Construct iterators from which to build the sub-string.
          // Note: The +1 while initializing beg is to skip past the "=" which follows the label name
          std::string::iterator
            beg = line.begin() + label_name.size() + 1 + label_index,
            end = (comma_index == std::string::npos) ? line.end() : line.begin() + comma_index;

          return std::string(beg, end);
        }

      // The label index was not found, return the empty string
      return std::string("");
    }



protected:
    std::ifstream m_file;
    Domain* m_domain;

};

} // namespace Oceane

#endif // MESHIO_H
