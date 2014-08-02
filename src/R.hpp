#ifndef R_HPP_INCLUDED
#define R_HPP_INCLUDED

#include <array>
#include <cassert>
#include <boost/format.hpp>


const int dims = 3;

/**
 * This is a vector representing the number of each nucleotide base in the current branch.
 */
class R {

private:
    /**
     * The total number of lineages.
     */
    int n;

    typedef std::array<int,dims> ValueArray;

    /**
     * The number of lineages with each nucleotide.
     */
    ValueArray values;

    int sum;



public:
    /**
     * Creates an R vector with a certian number of lineages and lineage counts.
     * @param n The number of lineages.
     * @param values The lineage counts.
     */
    R(int n, const std::array<int,dims>& values) : n(n), values(values){
        assert(isValidValues(n,values));

        sum = calculateSum(values);

    }


    static int calculateSum(const std::array<int,dims>& arr)
    {
        int sum = 0;
        for (int elem : arr)
            sum += elem;

        return sum;
    }

    /**
     * Gets the size of a matrix used to store all R for 1<=n<=m.
     * @param m The maximum number of lineages.
     * @return The size of the matrix.
     */
    static int getMatrixSize(int m){

        return calculateIndexedDimension(m,dims+1) -1;
    }

    /**
     * Gets the size of a matrix used to store all R with n total lineages.
     * @param n The total number of lineages.
     * @return The size of the matrix.
     */
    static int getMatrixSizeWithoutN(int n) {

        return calculateIndexedDimension(n,dims);
    }



    std::string toString() 
    {
        std::string resultString = str(boost::format("R{sum = %1%:")%sum);

        for (int type = 0; type < getNumberOfTypes(); type++)
        {
            resultString += (boost::format("%1% = %2%,") % getTypeName(type) % getNum(type)).str();
        }

        return  resultString + '}';
    }


    /**
     * Gets the "type name" of a type. These names are arbitrary and are mainly just used for debugging.
     * @param type The type
     * @return The name of that type.
     */
    std::string getTypeName(int type){
        switch (type)
        {
            case 0:
                return "r";
            case 1:
                return "g";
            case 2:
                return "y";
            case 3:
                return "t";

        }

        assert(false && "No such type");
    }


    /**
     * This calculates the number of indices for a vector of a given size and dimension.
     * How this works is not that clear.
     * @param index The size.
     * @param dimension The number of dimensions.
     * @return The number of indices.
     */
    static int calculateIndexedDimension(int index, int dimension)
    {
        int result = 1;

        for (int i = 1;i <= dimension; i++)
        {
            result  *= (index + i);
        }

        for (int i = 1; i <= dimension; i++) {

            result /=i;
        }

        return result;

    }

    /**
     * Gets the index of this R in a matrix with only one n value.
     * @return The index.
     */
    int getIndexWithoutN(){
        int total = 0;

        int totalCount = 0;
        for (int i = 0; i < dims; i++) {
            totalCount += values[i];
            total += calculateIndexedDimension(totalCount-1,i+1);
        }
        return total;
    }

    /**
     * Gets the index of this R in a matrix storing R's of varies n values.
     * @return The index.
     */
    int getIndex(){
        int total = calculateIndexedDimension(n-1,dims+1) -1;
        total += getIndexWithoutN();
        return total;
    }

    /**
     * Attempts to subtract one R from another.
     * Returns an optional as it can fail if it would result in negative size elements.
     * @param other The other R.
     * @return this - other
     */
    R subtract(const R& other)
    {
        int newN = n - other.n;
        std::array<int,dims> newValues = subtractValues(values,other.values);
        if (isValidValues(newN,newValues))
           return R(newN,newValues);
        else
            throw "Invalid subtractions";
    }

    /**
     * Checks if the values are valid for a given n.
     * @param n The number of lineages.
     * @param values The lineage counts.
     * @return True if they are valid.
     */
    static bool isValidValues(int n, const ValueArray& values)
    {
        if (n <=0)
            return false;

        int remaining = n;

        for (int a : values)
        {
            if (a<0)
                return false;
            remaining-=a;
        }

        return remaining >=0;
    }

    /**
     * Calculates vector a-b.
     * @return A new array holding a-b.
     */
    static ValueArray subtractValues(const ValueArray& a, const ValueArray& b)
    {
        ValueArray result;
        for (unsigned int i = 0; i < a.size(); i++) {
            result[i] = a[i] - b[i];
        }
        return result;
    }

    /*
    *
     * Returns the probability of selecting toSelect from this R.
     * Follows the multivariate hypergeometric distribution.
     * @param toSelect The R to select from.
     * @return The probablity of this selection.
     
    double getProbabilityOfSelecting(const R& toSelect)
    {
        double product = 1;
        for (int i =0; i < getNumberOfTypes(); i++)
        {
            product *= CombinatoricsUtils.binomialCoefficientDouble(getNum(i),toSelect.getNum(i));
        }
        product/= CombinatoricsUtils.binomialCoefficientDouble(n,toSelect.n);

        return product;
    }

    *
     * Calculates a factor for selection probabilities.
     * @return n!/(v[0]! * v[1]! .. * v[numberOfTypes]!).
     
    double getLikelihoodProduct()
    {
        double numerator = CombinatoricsUtils.factorialDouble(n);

        for (int i = 0; i < getNumberOfTypes(); i++) {
            numerator/= CombinatoricsUtils.factorialDouble(getNum(i));
        }

        if (Double.isInfinite(numerator) || Double.isNaN(numerator))
            throw new RuntimeException("GetLikelihoodproduct is not returning a finite value.");

        return numerator;
    }
    */


    /**
     * Gets the index for R in a square matrix of a size.
     * @param size The size of the square matrix.
     * @return The index.
     */
    int getBoxIndex(int size)
    {
        int sum =n;

        for (int i = 0; i < dims; i++)
        {
            sum *= size;
            sum += getNum(i);
        }
        return sum;
    }

    /**
     * Returns an iterable over all values of R for a given n.
     * @param n The number of lineages.
     * @return An iterable to iterate over.
     
    static Iterable<R> loopOver(final int n) {
        return new Iterable<R>()
        {
            @Override
            public Iterator<R> iterator()
            {
                return new Iterator<R>()
                {

                    R last = null;

                    @Override
                    public boolean hasNext()
                    {
                        return last == null || last.values[dims - 1] != n;
                    }

                    @Override
                    public R next()
                    {
                        if (last == null)
                            last = new R(n, new int[dims]);
                        else
                        {
                            tryAdvance(0);
                        }
                        last.sum = calculateSum(last.values);
                        return last;
                    }

                    @Override
                    public void remove()
                    {
                        throw new UnsupportedOperationException();
                    }

                    private void tryAdvance(int i)
                    {
                        if (i == -1)
                            throw new RuntimeException("Iterator went over limit");


                        if (last.sum == n)
                        {
                            last.values[i] = 0;
                            last.sum = calculateSum(last.values);
                            tryAdvance(i + 1);
                        } else
                            last.values[i]++;

                    }
                };
            }
        };

    }
    */

    /**
     * @return The number of types for this R.
     */
    static int getNumberOfTypes()
    {
        return dims+1;
    }


    /**
     * Gets the number of lineages of a certain type from the R vector.
     * @param type The type
     * @return The number of lineages of that type.
     */
    int getNum(int type)
    {
        if (type != dims)
            return values[type];

        else
            return n-sum;
    }


    /**
     * Assuming the given R value represents only one lineage, returns the type of that lineage.
     * @return The type
     */
    int getType() {
        if (n!= 1)
            throw "A type should only have one allele";

        for (int i = 0; i < getNumberOfTypes();i++)
        {
            if (getNum(i) == 1)
                return i;
        }

        throw "Should be either red or green";
    }

    /**
     * Equivalent to transition(fromType,toType).getIndex();
     */
    int transitionIndex(int fromType, int toType) {

        if (fromType != dims)
            values[fromType] -=1;

        if (toType != dims)
            values[toType] += 1;

        int result =  getIndex();

        if (fromType != dims)
            values[fromType] +=1;

        if (toType != dims)
            values[toType] -= 1;

        return result;
    }

    /**
     * Simulates a mutation from the fromType to the toType.
     * @param fromType The source lineage type.
     * @param toType The destination lineage type.
     * @return A new R representing after the mutation.
     */
    R transition(int fromType, int toType) {

        ValueArray newValues(values);

        if (fromType != dims)
            newValues[fromType] -=1;

        if (toType != dims)
            newValues[toType] += 1;

        return R(n,newValues);
    }

    /**
     * Simulates coalescing two of the lineages into one.
     * @param coalesceType The type to coalesce.
     * @return A new R after coalesing.
     */
    R coalesce(int coalesceType) {

        ValueArray newValues(values);

        if (coalesceType != dims)
            newValues[coalesceType] -=1;

        return R(n-1,newValues);
    }

    /**
     * Equivalent to split(splitType).getIndex()
     */
    int splitIndex(int splitType) {

        if (splitType != dims)
            values[splitType] += 1;

        n+= 1;

        int result = getIndex();
        if (splitType != dims)
            values[splitType] -= 1;

        n-=1;

        return result;
    }

    /**
     * Simulates a lineage splitting into two. (IE, the opposite of coalescing).
     * @param splitType The type to split.
     * @return A new R after this split.
     */
    R split(int splitType) {
        ValueArray newValues(values);

        if (splitType != dims)
            newValues[splitType] += 1;

        return R(n+1,newValues);
    }

    friend class RIterator;


           



    /*
    static void main()
    {

        for (int n = 1; n < 100; n++)
        {
            System.out.println(n+","+R.getMatrixSize(n));
        }


        Map<Integer,String> seenIndices = new HashMap<>();
        int count = 0;

        for (int n = 1; n <=20 ;n++) {


            for (R r : R.loopOver(n)) {
                count++;
                int firstIndex = r.getIndex();

                if (seenIndices.containsKey(firstIndex)) {
                    throw new RuntimeException("Fail");
                }

                seenIndices.put(firstIndex, r.toString());

                System.out.println(firstIndex + "," + r.getIndexWithoutN() + " , " + r);

            }
        }
        System.out.println(count);
        System.out.println(seenIndices.size());

    }*/

};

    class RIterator
    {
        private:
            int n;
            R last;
            bool invalid;
        public:
            RIterator(int n) : n(n), last(n,R::ValueArray()), invalid(false)
            {
            }

            RIterator() : last(1,R::ValueArray()),invalid(true)
            {
            }

            void tryAdvance(int i)
            {
                if (i == -1)
                {
                    invalid = true;
                    return;
                }

                if (last.sum == n)
                {
                    last.sum -= last.values[i];
                    last.values[i] = 0;
                    tryAdvance(i - 1);
                } 
                else
                {
                    last.values[i]++;
                    last.sum++;
                }

            }

            R& operator*()
            {
                assert(!invalid && "Iterator is not valid");
                return last;
            }

            bool operator != (const RIterator& other)
            {
                return !(*this == other);
            }

            bool operator == (const RIterator& other)
            {
                return invalid && other.invalid;
            }
                    

            RIterator& operator++()
            {
                assert(!invalid && "Iterator is not valid");
                tryAdvance(dims-1);
                return *this;
            }
     };

    class IterateUpTo
    {
        private:
            const int n;
        public:
            IterateUpTo(int n): n(n)
            {
            }

            RIterator begin()
            {
                return RIterator(n);
            }

            RIterator end()
            {
                return RIterator();
            }
    };

    IterateUpTo loopOver(int n)
    {
        return IterateUpTo(n);
    }

#endif
