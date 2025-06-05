#ifndef BASE_LEARNED_MODEL
#define BASE_LEARNED_MODEL

#include <vector>

template <typename KeyType>
class baseLearnedModel
{
    public:
        virtual void buildIndex(KeyType data[], long startPos, long endPos) = 0;
        virtual long getPosition(KeyType key) = 0;
};
#endif