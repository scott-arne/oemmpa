#include "oemmpa/Analyzer.h"
#include "oemmpa/MatchedPair.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

/// \brief The seven-field natural key used by the golden pair dumps.
std::string natural_key(const OEMMPA::MatchedPair& pair) {
    std::ostringstream key;
    key << pair.GetSourceExternalId() << '\t'
        << pair.GetTargetExternalId() << '\t'
        << pair.GetConstantSmiles() << '\t'
        << pair.GetSourceVariableSmiles() << '\t'
        << pair.GetTargetVariableSmiles() << '\t'
        << pair.GetCutCount() << '\t'
        << pair.GetHeavyAtomDelta() << '\t'
        << pair.GetHeavyBondDelta();
    return key.str();
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: " << argv[0] << " <input.smi> <output.tsv>\n";
        return 2;
    }
    std::ifstream input(argv[1]);
    if (!input) {
        std::cerr << "cannot read " << argv[1] << "\n";
        return 2;
    }

    OEMMPA::Analyzer analyzer("wizepairz");
    std::string smiles;
    std::string identifier;
    while (input >> smiles >> identifier) {
        analyzer.AddMolecule(smiles, identifier);
    }
    analyzer.Analyze();

    std::vector<std::string> keys;
    for (const OEMMPA::MatchedPair& pair : analyzer.GetPairs()) {
        keys.push_back(natural_key(pair));
    }
    // Sorted so the comparison is order-independent, but deliberately NOT
    // deduplicated: a pair emitted twice on one side must not hide behind an
    // identical key on the other.
    std::sort(keys.begin(), keys.end());

    std::ofstream output(argv[2]);
    for (const std::string& key : keys) {
        output << key << '\n';
    }
    if (!output) {
        std::cerr << "cannot write " << argv[2] << "\n";
        return 1;
    }
    return 0;
}
