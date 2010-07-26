#ifndef _CORE_FACTORY_HH
#define _CORE_FACTORY_HH
#include <map>

/**
 * Object factory
 * based on Alexandrescu, Andrei "Modern C++ Design"
 */
template<class BaseClass, typename CreationFunction, typename Identifier>
class Factory
{
  typedef std::map<Identifier, CreationFunction> Registry;
public:
  Factory() {}

  /**
   * register a class
   * @param id identifier for this class
   * @param c  function that creates a new instance of this class
   * @return class registered
   */
  bool registerClass(const Identifier &id, CreationFunction c) {
    if(registry_.find(id) == registry_.end()) {
      registry_.insert(typename Registry::value_type(id, c));
      return true;
    }
    else
      return false;
  }

  /**
   * Create an instance of the class which is identified by @c id
   * If no such classes exists 0 is returned !
   */
  BaseClass* getObject(const Identifier &id) const{
    CreationFunction create = getCreationFunction(id);
    if(create)
      return create();
    else
      return 0;
  }

  /**
   * Create an instance of the class which is identified by @c id
   * If no such classes exists 0 is returned !
   */
  template<class Param>
  BaseClass* getObject(const Identifier &id, const Param& p) {
    CreationFunction create = getCreationFunction(id);
    if(create)
      return create(p);
    else
      return 0;
  }

  /**
   * Create an instance of the class which is identified by @c id
   * If no such classes exists 0 is returned !
   */
  template<class Param1, class Param2>
  BaseClass* getObject(const Identifier &id, const Param1& p1, const Param2& p2) {
    CreationFunction create = getCreationFunction(id);
    if(create)
      return create(p1,p2);
    else
      return 0;
  }



protected:
  CreationFunction getCreationFunction(const Identifier &id) const {
    typename Registry::const_iterator item = registry_.find(id);
    if(item != registry_.end())
      return item->second;
    else
      return 0;
  }

  Registry registry_;
};


#endif
