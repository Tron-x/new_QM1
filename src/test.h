#include "eft_calculator.hpp"

using namespace std;

vector<double> test_inp(const string& filename, qm::EFTCalculator& m_calculator) {
    ifstream ifs(filename);
    if (!ifs.good()) {
        cerr << "ERROR: Opening file '" << filename << "' failed.\n";
        exit(EXIT_FAILURE);
    }

    vector<string> lines;
    string line;
    while (!ifs.eof()) {
        getline(ifs, line);
        common::trim(line);

        if (line == "" || line == "$END" || line == "$end") {
            break;
        }
        lines.push_back(line);
    }

    vector<vector<double>> coors;
    for (int idx = lines.size() - 6; idx < lines.size(); idx++) {
        auto fields = common::tokenizer(lines[idx], ' ');
        assert(fields.size() == 5);
        coors.push_back(
                {stod(fields[2]), stod(fields[3]), stod(fields[4])});
    }
    ifs.close();

    auto eft = m_calculator.eval(coors);
    return eft;
}