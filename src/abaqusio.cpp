#include "abaqusio.h"
#include <memory>

#include "elementfactory.h"
namespace Oceane {


AbaqusIO::AbaqusIO(std::string szFile,Domain* domain) :m_domain(domain) {
    std::string lastThree = szFile.substr(szFile.length() - 3);
    if(lastThree == "inp") {
        m_file.open(szFile,std::ios::in);
        std::cout<<szFile<<"\tOpened Successfully\n";
        if(!m_file)
            std::cout<<"Coudn't open the file\n";
        this->read();
    }
    if (lastThree == "csv") {
        prefix_file_name = szFile.substr(0, szFile.length() - 4);
        this->read_csv();
    }
}

void AbaqusIO::read_csv()
{
    // NodeMap m_nodemap;
    std::ifstream inpFile(prefix_file_name + "_node_data.csv");
    if (inpFile.is_open()) {
        std::cout << "File opened successfully" << std::endl;
    }

    Index nodeId;
    double x, y;

    std::string nodeId_s, x_s, y_s;

    std::string first_line;
    getline(inpFile, first_line); // skip the first line
    while (getline(inpFile, first_line)) {
        getline(inpFile, nodeId_s, ',');
        getline(inpFile, x_s, ',');
        getline(inpFile, y_s);

        nodeId = (Index)stoi(nodeId_s);
        x = stod(x_s);
        y = stod(y_s);

        m_nodemap[nodeId] = std::vector<double>{x,y};
    }

    // ElementMap m_elementmap; only 2 nodes per member 
    std::ifstream inpFile(prefix_file_name + "_member_data.csv");
    if (inpFile.is_open()) {
        std::cout << "File opened successfully" << std::endl;
    }

    Index nodeId, m1, m2;
    double area, matparam1, matparam2, matparam3, matparam4;
    std::map<Index,std::vector<Index> > elements;

    std::string material, nodeId_s, m1_s, m2_s, area_s, matparam1_s, matparam2_s, matparam3_s, matparam4_s;

    std::string first_line;
    getline(inpFile, first_line); // skip the first line
    while (getline(inpFile, first_line)) {
        getline(inpFile, material, ',');
        getline(inpFile, nodeId_s, ',');
        getline(inpFile, m1_s, ',');
        getline(inpFile, m2_s, ',');
        getline(inpFile, area_s, ',');
        getline(inpFile, matparam1_s, ',');
        getline(inpFile, matparam2_s, ',');
        getline(inpFile, matparam3_s, ',');
        getline(inpFile, matparam4_s);

        nodeId = (Index)stoi(nodeId_s);
        m1 = (Index) stoi(m1_s);
        m2 = (Index) stoi(m2_s);

        elements[nodeId] = std::vector<Index> {m1, m2};
        m_elementmap["ELEMENT"] = elements;
    }

    // material = "model" as per the example
    area = stod(area_s);
    matparam1 = stod(matparam1_s);
    matparam2 = stod(matparam2_s);
    matparam3 = stod(matparam3_s);
    matparam4 = stod(matparam4_s);

    // LoadCsv   m_loads_csv; Might need to change the naming convention
    std::ifstream inpFile(prefix_file_name + "_force_data.csv");
    if (inpFile.is_open()) {
        std::cout << "File opened successfully" << std::endl;
    }

    Index forceId;
    double x, y;

    std::string forceId_s, Fx, Fy;

    std::string first_line;
    getline(inpFile, first_line); // skip the first line
    while (getline(inpFile, first_line)) {
        getline(inpFile, forceId_s, ',');
        getline(inpFile, Fx, ',');
        getline(inpFile, Fy);

        forceId = (Index) stoi(forceId_s);
        x = stod(Fx);
        y = stod(Fy);

        m_loads_csv[forceId] = std::vector<double>{x,y};
    }

    // FixityCsv m_fixity_csv;
    std::ifstream inpFile(prefix_file_name + "_support.csv");
    if (inpFile.is_open()) {
        std::cout << "File opened successfully" << std::endl;
    }

    Index fixityId;
    double x, y;

    std::string fixityId_s, ux, uy;

    std::string first_line;
    getline(inpFile, first_line); // skip the first line
    while (getline(inpFile, first_line)) {
        getline(inpFile, fixityId_s, ',');
        getline(inpFile, ux, ',');
        getline(inpFile, uy);

        fixityId = (Index) stoi(fixityId_s);
        x = stod(ux);
        y = stod(uy);

        m_fixity_csv[fixityId] = std::vector<double>{x,y};
    }

    // BoundCondn m_bcmap;
    std::ifstream inpFile(prefix_file_name + "_boundary_conditions.csv");
    if (inpFile.is_open()) {
        std::cout << "File opened successfully" << std::endl;
    }

    Index bcId;
    double x, y;

    std::string bcId_s, ux, uy;

    std::string first_line;
    getline(inpFile, first_line); // skip the first line
    while (getline(inpFile, first_line)) {
        getline(inpFile, bcId_s, ',');
        getline(inpFile, ux, ',');
        getline(inpFile, uy);

        bcId = (Index) stoi(bcId_s);
        x = stod(ux);
        y = stod(uy);

        m_bcmap[bcId] = std::vector<double>{x,y};
    }

}

void AbaqusIO::read()
{
        //read the file line by line
        std::string s;
        while(std::getline(m_file,s))
        {
            //convert all letters into upper.
            std::string upper(s);
            std::transform(upper.begin(),upper.end(),upper.begin(),toupper);
            upper.erase(std::remove_if(upper.begin(), upper.end(), isspace), upper.end());
            if(upper.find("*NODE") != std::string::npos)
            {
                readNodes();
            }

            if(upper.find("*ELEMENT") != std::string::npos)
            {
                readElements(upper);
            }

            if(upper.find("*NSET") != std::string::npos)
            {
                readNset(upper);

            }

            if(upper.find("*ELSET") != std::string::npos)
            {
                readElset(upper);

            }
            if(upper.find("*BOUNDARY")!= std::string::npos )
            {
                readFixity();
            }

            if(upper.find("*DLOAD") != std::string::npos)
            {
                readLoad();
            }

        }

        this->addToDomain();
        m_domain->Init();
        std::cout<<"Domain Initiated\n";
}

void AbaqusIO:: readElements(std::string upper)
{

    std::string  elType = parse_label(upper,"TYPE");
    Index id; char c;
    //below loop for quad4 element

    std::map<Index,std::vector<Index> > elements;

    while(m_file.peek() != '*' and m_file.peek() != EOF)
    {

        if(m_file.peek() == '\n')
        {
            std::string dummy;
            std::getline(m_file,dummy); continue;    // skip the blank line.
        }

        // TODO : ENABLE VECTOR TO STORE IDS, IRRESPECTIVE OF SIZE.
        Index n1,n2,n3,n4;
        if(m_file>> id >> c >> n1>> c >> n2>> c >> n3>> c >>n4)
        {
            std::vector<Index> nodeIds= {n1,n2,n3,n4};
            elements[id]= nodeIds;

            std::string dummy;
            std::getline(m_file,dummy);
            continue;  // dummy string to collect junk data is not used
        }
    }

    m_elementmap[elType] =elements;
}



void AbaqusIO:: readNodes()
{
    //TODO: Make dimension independent. add z coordinate.

    std::string line;
    char c;                         //temp.variable for recieving ','.
    Oceane::Index id;
    double x,y;
    std::vector<double> coords;

    //The program will fail if the input file contains blank lines.
    while(m_file.peek() != '*' and m_file.peek() != EOF )
    {

        if(m_file.peek() == '\n')              //if the line is blank skip to next line.
        {
            std::string skipline;
            std::getline(m_file,skipline); continue;
        }

        std::getline(m_file,line);
        line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
        std::stringstream ss(line);

        ss>>id>>c>>x>>c>>y;

        //update nodeId - coordinate map i.e m_nodemap.
        m_nodemap[id] = std::vector<double>{x,y};
    }
}



void AbaqusIO::readElset(std::string line)
{
    std::string elset_name = parse_label(line,"ELSET");
    while(m_file.peek() != '*' and m_file.peek() != EOF )
    {
        if(m_file.peek() == '\n')
        {
            std::string dummy;
            std::getline(m_file,dummy);continue;        // skip the blank line.
        }
        std::string ss; Index id;
        std::getline(m_file,ss);
        std::transform(ss.begin(),ss.end(),ss.begin(),toupper);
        std::stringstream str(ss);
        std::vector<Index> element_ids{};
        while(str>>id)
        {
            element_ids.push_back(id);
            if(str.peek() == ',' or str.peek() == ' ')
            {
                str.ignore();
            }

        }
        m_elsetmap[elset_name] = element_ids;
    }
}



void AbaqusIO::readNset(std::string line)
{
    std::string nodeset_name = parse_label(line,"NSET");
    std::vector<Index> id_s;
    while(m_file.peek() != '*' and m_file.peek() != EOF )
    {
        if(m_file.peek() == '\n')
        {
            std::string dummy;
            std::getline(m_file,dummy);continue;        // skip the blank line.
        }
        Index id;
        std::string ss;
        std::getline(m_file,ss);
        std::stringstream str(ss);

        while(str>>id)
        {
            id_s.push_back(id);
            if(str.peek() == ',' or str.peek() == ' ')
            {
                str.ignore();
            }
        }
    }
    m_nsetmap[nodeset_name] =id_s;
}


void AbaqusIO::addToDomain()
{
//    adding nodes To domain
    for(auto& iter:m_nodemap)
    {
        m_domain->add_Nodeptr(std::make_shared<Node>(iter.first,iter.second));
    }

//    adding elements to domain
    ElementFactory factory;
    for(auto& iter:m_elementmap)
    {
        std::string ElementType = iter.first;
        for(auto& elem:iter.second)
        {
            Index elemId =elem.first;
            auto nodes = m_domain->getNodes_at(elem.second);
            m_domain->add_Elementptr(factory.createElement(ElementType ,elemId,nodes));
        }
    }

//    adding Nset map
    std::map<std::string,std::vector<Nodeptr> > nsetmap;
    for(auto& iter:m_nsetmap)
    {
        std::vector<Nodeptr> _nodes{};
        for(auto& id:iter.second)
        {
            _nodes.push_back(m_domain->getNodeptr(id));

        }
        nsetmap[iter.first] = _nodes;
    }
    m_domain->setNsetmap(nsetmap);

//    adding Elset map
    std::map<std::string, std::vector<Elementptr> > elsetmap;
    for(auto& iter:m_elsetmap)
    {
        std::vector<Elementptr> _elements{};
        for(auto& id:iter.second)
        {
            _elements.push_back(m_domain->getElementptr(id));
        }
        elsetmap[iter.first]=_elements;
    }
    m_domain->setElsetmap(elsetmap);

//    adding load
    m_domain->setLoads(m_loads);

//    adding diplacement constraint
    m_domain->setFixity(m_fixity);

}




/* read lines  ELSET, <faceNo>, <Load>  */
void AbaqusIO::readLoad()
{
    while(m_file.peek() != '*' and m_file.peek() != EOF )
    {
        if(m_file.peek() == '\n')
        {
            std::string dummy;
            std::getline(m_file,dummy);continue;        // skip the blank line.
        }

        std::string ss;   // line containg elset,face,pressure
        std::getline(m_file,ss);

        // Remove whitespace.
        strip_ws(ss);
        std::transform(ss.begin(),ss.end(),ss.begin(),toupper);

        std::string elset_name, szface; double pressure;
        {
            char comma =','; char c;
            std::stringstream ss_stream(ss);
            {
                int flag =0;
                while(ss_stream>>c)
                {

                    if(c==comma)
                       {
                           flag+=1;
                           if(flag ==2)
                               break;
                    }

                }
                ss_stream>>pressure;        //get pressure load.
            }

            size_t first_comma_index, second_comma_index;
             first_comma_index = ss.find(",");
            second_comma_index = ss.find(",",first_comma_index);
            {
                std::string::iterator
                        beg =ss.begin(),
                        end =ss.begin()+first_comma_index;
                elset_name= std::string(beg,end);
            }
            {
               std::string::iterator
                        beg =ss.begin() +first_comma_index+1,
                        end =ss.begin() +first_comma_index+3;;
                szface =std::string(beg,end);
            }   //get element set and load orientation

        }

        unsigned int nface=0;
        if(szface=="P2")
            nface =1;
        else if(szface=="P3")
            nface =2;
        else if(szface =="P4")
            nface =3;

        Oceane::Load _load{elset_name,nface,pressure};
        m_loads.push_back(_load);
    }
}



/* read lines  ELSET, <faceNo> */
//TODO: currently constraints in all direction is put in both nodes of the face
//     restrictions can be made more generalized.
void AbaqusIO::readFixity()
{
    while(m_file.peek() != '*' and m_file.peek() != EOF )
    {
        if(m_file.peek() == '\n')
        {
            std::string dummy;
            std::getline(m_file,dummy);continue;        // skip the blank line.
        }

        std::string ss;   // line containg elset,face,pressure
        std::getline(m_file,ss);

        // Remove whitespace.
        strip_ws(ss);
        std::transform(ss.begin(),ss.end(),ss.begin(),toupper);

        std::string elset_name, szface; double pressure;
        {

            size_t first_comma_index, second_comma_index;
             first_comma_index = ss.find(",");
            second_comma_index = ss.find(",",first_comma_index);
            {
                std::string::iterator
                        beg =ss.begin(),
                        end =ss.begin()+first_comma_index;
                elset_name= std::string(beg,end);
            }
            {
               std::string::iterator
                        beg =ss.begin() +first_comma_index+1,
                        end =ss.begin() +first_comma_index+3;;
                szface =std::string(beg,end);
            }   //get element set and load orientation

        }

        unsigned int nface=0;
        if(szface=="P2")
            nface =1;
        else if(szface=="P3")
            nface =2;
        else if(szface =="P4")
            nface =3;

        Oceane::Fixity _fix{elset_name,nface};
        m_fixity.push_back(_fix);
    }
}


} // namespace Oceane
