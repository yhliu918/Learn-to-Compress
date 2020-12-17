/**
 * This code is released under the
 * Apache License Version 2.0 http://www.apache.org/licenses/.
 *
 * (c) Daniel Lemire, http://lemire.me/en/
 */

#ifndef CODECFACTORY_H_
#define CODECFACTORY_H_

#include "common.h"
#include "codecs.h"
#include "util.h"
#include "FOR.h"
#include "piecewise.h"
#include "piecewise_fix.h"
#include "maskvbyte.h"



namespace Codecset {

typedef std::map<std::string, IntegerCODEC*> CodecMap;

/**
 * This class is a convenience class to generate codecs quickly.
 * It cannot be used safely in a multithreaded context where
 * each thread should have a different IntegerCODEC.
 */
class CODECFactory {
public:
  static CodecMap scodecmap;

  // hacked for convenience
  static std::vector<IntegerCODEC*> allSchemes() {
    std::vector<IntegerCODEC*> ans;
    for (auto i = scodecmap.begin(); i != scodecmap.end(); ++i) {
      ans.push_back(i->second);
    }
    return ans;
  }

  static std::vector<std::string> allNames() {
    std::vector<std::string> ans;
    for (auto i = scodecmap.begin(); i != scodecmap.end(); ++i) {
      ans.push_back(i->first);
    }
    return ans;
  }

  static IntegerCODEC* &getFromName(std::string name) {
    if (scodecmap.find(name) == scodecmap.end()) {
      std::cerr << "name " << name << " does not refer to a CODEC."
                << std::endl;
      std::cerr << "possible choices:" << std::endl;
      for (auto i = scodecmap.begin(); i != scodecmap.end(); ++i) {
        std::cerr << static_cast<std::string>(i->first)
                  << std::endl; // useless cast, but just to be clear
      }
      std::cerr << "for now, I'm going to just return 'copy'" << std::endl;
      return scodecmap["copy"];
    }
    return scodecmap[name];
  }
};

// C++11 allows better than this, but neither Microsoft nor Intel support C++11
// fully.
static inline CodecMap initializefactory() {
  CodecMap map;

  map["FOR"]= new Codecset::FOR();
  
  map["piecewise"]= new piecewise();
  map["piecewise_fix"]= new piecewise_fix();
  //map["BP32"] = new CompositeCodec<BP32, VariableByte>;
  map["MaskVByte"] = new Codecset::MaskVByte();
  map["copy"] = new Codecset::JustCopy();
  return map;
}

CodecMap CODECFactory::scodecmap = initializefactory();

} // namespace FastPFor

#endif /* CODECFACTORY_H_ */
