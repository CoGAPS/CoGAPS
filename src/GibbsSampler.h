#ifndef __COGAPS_GIBBS_SAMPLER_H__
#define __COGAPS_GIBBS_SAMPLER_H__

#include "Archive.h"
#include "AtomicDomain.h"
#include "GapsAssert.h"
#include "ProposalQueue.h"
#include "data_structures/Matrix.h"
#include "math/Algorithms.h"
#include "math/Random.h"

#include <algorithm>

// forward declarations needed for friend classes/functions
class AmplitudeGibbsSampler;
class PatternGibbsSampler;
class GapsStatistics;

template <class T, class MatA, class MatB>
class GibbsSampler;

template <class T, class MatA, class MatB>
inline Archive& operator<<(Archive &ar, GibbsSampler<T, MatA, MatB> &samp);

template <class T, class MatA, class MatB>
inline Archive& operator>>(Archive &ar, GibbsSampler<T, MatA, MatB> &samp);

/*************************** GIBBS SAMPLER INTERFACE **************************/

template <class T, class MatA, class MatB>
class GibbsSampler
{
private:

    friend class GapsStatistics;

protected:

    MatB mDMatrix;
    MatB mSMatrix;
    MatB mAPMatrix;

    MatA mMatrix;
    MatB* mOtherMatrix;

    ProposalQueue mQueue;
    AtomicDomain mDomain;

    float mLambda;
    float mMaxGibbsMass;
    float mAnnealingTemp;
    
    unsigned mNumRows;
    unsigned mNumCols;
    uint64_t mBinSize;

    float mAvgQueue;
    float mNumQueues;

    T* impl();

    void processProposal(AtomicProposal *prop);

    void birth(AtomicProposal *prop);
    void death(AtomicProposal *prop);
    void move(AtomicProposal *prop);
    void exchange(AtomicProposal *prop);

    bool updateAtomMass(Atom *atom, float delta);
    void acceptExchange(AtomicProposal *prop, float d1, unsigned r1,
        unsigned c1, unsigned r2, unsigned c2);

    std::pair<float, bool> gibbsMass(AlphaParameters alpha, GapsRng *rng);
    std::pair<float, bool> gibbsMass(AlphaParameters alpha, float m1, float m2,
        GapsRng *rng);

public:

    template <class DataType>
    GibbsSampler(const DataType &data, bool transpose, unsigned nPatterns,
        bool amp, bool partitionRows, const std::vector<unsigned> &indices);

    template <class DataType>
    void setUncertainty(const DataType &unc, bool transposeData,
        bool partitionRows, const std::vector<unsigned> &indices);
    
    void setSparsity(float alpha, bool singleCell);
    void setMaxGibbsMass(float max);
    void setAnnealingTemp(float temp);

    void setMatrix(const Matrix &mat);

    void update(unsigned nSteps, unsigned nCores);

    unsigned dataRows() const;
    unsigned dataCols() const;

    float chi2() const;
    uint64_t nAtoms() const;

    #ifdef GAPS_DEBUG
    float getAvgQueue() const;
    bool internallyConsistent();
    #endif

    // serialization
    friend Archive& operator<< <T, MatA, MatB> (Archive &ar, GibbsSampler &sampler);
    friend Archive& operator>> <T, MatA, MatB> (Archive &ar, GibbsSampler &sampler);
};

class AmplitudeGibbsSampler : public GibbsSampler<AmplitudeGibbsSampler, ColMatrix, RowMatrix>
{
private:

    friend class GibbsSampler;
    friend class PatternGibbsSampler;

    unsigned getRow(uint64_t pos) const;
    unsigned getCol(uint64_t pos) const;
    bool canUseGibbs(unsigned row, unsigned col) const;
    bool canUseGibbs(unsigned r1, unsigned c1, unsigned r2, unsigned c2) const;
    void updateAPMatrix(unsigned row, unsigned col, float delta);

    AlphaParameters alphaParameters(unsigned row, unsigned col);
    AlphaParameters alphaParameters(unsigned r1, unsigned c1, unsigned r2,
        unsigned c2);

    float computeDeltaLL(unsigned row, unsigned col, float mass);
    float computeDeltaLL(unsigned r1, unsigned c1, float m1, unsigned r2,
        unsigned c2, float m2);

public:

    template <class DataType>
    AmplitudeGibbsSampler(const DataType &data, bool transposeData,
        unsigned nPatterns, bool partitionRows,
        const std::vector<unsigned> &indices);

    void sync(PatternGibbsSampler &sampler);
    void recalculateAPMatrix();
};

class PatternGibbsSampler : public GibbsSampler<PatternGibbsSampler, RowMatrix, ColMatrix>
{
private:

    friend class GibbsSampler;
    friend class AmplitudeGibbsSampler;

    unsigned getRow(uint64_t pos) const;
    unsigned getCol(uint64_t pos) const;
    bool canUseGibbs(unsigned row, unsigned col) const;
    bool canUseGibbs(unsigned r1, unsigned c1, unsigned r2, unsigned c2) const;
    void updateAPMatrix(unsigned row, unsigned col, float delta);

    AlphaParameters alphaParameters(unsigned row, unsigned col);
    AlphaParameters alphaParameters(unsigned r1, unsigned c1, unsigned r2,
        unsigned c2);

    float computeDeltaLL(unsigned row, unsigned col, float mass);
    float computeDeltaLL(unsigned r1, unsigned c1, float m1, unsigned r2,
        unsigned c2, float m2);

public:

    template <class DataType>
    PatternGibbsSampler(const DataType &data, bool transposeData,
        unsigned nPatterns, bool partitionRows,
        const std::vector<unsigned> &indices);

    void sync(AmplitudeGibbsSampler &sampler);
    void recalculateAPMatrix();
};

/******************** IMPLEMENTATION OF TEMPLATED FUNCTIONS *******************/

template <class DataType>
AmplitudeGibbsSampler::AmplitudeGibbsSampler(const DataType &data,
bool transposeData, unsigned nPatterns, bool partitionRows,
const std::vector<unsigned> &indices)
    :
GibbsSampler(data, transposeData, nPatterns, true, partitionRows, indices)
{}

template <class DataType>
PatternGibbsSampler::PatternGibbsSampler(const DataType &data,
bool transposeData, unsigned nPatterns, bool partitionRows,
const std::vector<unsigned> &indices)
    :
GibbsSampler(data, transposeData, nPatterns, false, partitionRows, indices)
{}

template <class T, class MatA, class MatB>
template <class DataType>
GibbsSampler<T, MatA, MatB>::GibbsSampler(const DataType &data,
bool transposeData, unsigned nPatterns, bool amp, bool partitionRows,
const std::vector<unsigned> &indices)
    :
mDMatrix(data, transposeData, partitionRows, indices),
mSMatrix(mDMatrix.pmax(0.1f, 0.1f)),
mAPMatrix(mDMatrix.nRow(), mDMatrix.nCol()),
mMatrix(amp ? mDMatrix.nRow() : nPatterns, amp ? nPatterns : mDMatrix.nCol()),
mOtherMatrix(NULL), mQueue(amp ? mDMatrix.nRow() : mDMatrix.nCol(), nPatterns),
mDomain(mMatrix.nRow() * mMatrix.nCol()), mLambda(0.f),
mMaxGibbsMass(100.f), mAnnealingTemp(1.f), mNumRows(mMatrix.nRow()),
mNumCols(mMatrix.nCol()), mAvgQueue(0.f), mNumQueues(0.f)
{
    // calculate atomic domain size
    mBinSize = std::numeric_limits<uint64_t>::max()
        / static_cast<uint64_t>(mNumRows * mNumCols);

    // default sparsity parameters
    setSparsity(0.01, false);
}

template <class T, class MatA, class MatB>
template <class DataType>
void GibbsSampler<T, MatA, MatB>::setUncertainty(const DataType &unc,
bool transpose, bool partitionRows, const std::vector<unsigned> &indices)
{
    mSMatrix = MatB(unc, transpose, partitionRows, indices);
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::setSparsity(float alpha, bool singleCell)
{
    mQueue.setAlpha(alpha);

    float meanD = singleCell ? gaps::algo::nonZeroMean(mDMatrix) :
        gaps::algo::mean(mDMatrix);

    unsigned nPatterns = mDMatrix.nRow() == mMatrix.nRow() ? mMatrix.nCol() :
        mMatrix.nRow();

    mLambda = alpha * std::sqrt(nPatterns / meanD);
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::setMaxGibbsMass(float max)
{
    mMaxGibbsMass = max;
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::setAnnealingTemp(float temp)
{
    mAnnealingTemp = temp;
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::setMatrix(const Matrix &mat)
{   
    mMatrix = mat;
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::update(unsigned nSteps, unsigned nCores)
{
    unsigned n = 0;
    while (n < nSteps)
    {
        // populate queue, prepare domain for this queue
        //mQueue.populate(mDomain, nSteps - n);
        mQueue.populate(mDomain, 1);
        mDomain.resetCache(mQueue.size());
        n += mQueue.size();
        
        // update average queue count
        #ifdef GAPS_DEBUG
        mNumQueues += 1.f;
        mAvgQueue *= (mNumQueues - 1.f) / mNumQueues;
        mAvgQueue += mQueue.size() / mNumQueues;
        #endif

        // process all proposed updates
        #pragma omp parallel for num_threads(nCores)
        for (unsigned i = 0; i < mQueue.size(); ++i)
        {
            processProposal(&mQueue[i]);
        }

        mDomain.flushCache();
        mQueue.clear();
    }
}

template <class T, class MatA, class MatB>
unsigned GibbsSampler<T, MatA, MatB>::dataRows() const
{
    return mDMatrix.nRow();
}

template <class T, class MatA, class MatB>
unsigned GibbsSampler<T, MatA, MatB>::dataCols() const
{
    return mDMatrix.nCol();
}

template <class T, class MatA, class MatB>
float GibbsSampler<T, MatA, MatB>::chi2() const
{   
    return 2.f * gaps::algo::loglikelihood(mDMatrix, mSMatrix, mAPMatrix);
}
  
template <class T, class MatA, class MatB>
uint64_t GibbsSampler<T, MatA, MatB>::nAtoms() const
{   
    return mDomain.size();
}

#ifdef GAPS_DEBUG
template <class T, class MatA, class MatB>
float GibbsSampler<T, MatA, MatB>::getAvgQueue() const
{
    return mAvgQueue;
}

template <class T, class MatA, class MatB>
bool GibbsSampler<T, MatA, MatB>::internallyConsistent()
{
    return true;
}
#endif

template <class T, class MatA, class MatB>
Archive& operator<<(Archive &ar, GibbsSampler<T, MatA, MatB> &sampler)
{
    ar << sampler.mMatrix << sampler.mAPMatrix << sampler.mQueue <<
        sampler.mDomain << sampler.mLambda << sampler.mMaxGibbsMass <<
        sampler.mAnnealingTemp << sampler.mNumRows << sampler.mNumCols <<
        sampler.mBinSize << sampler.mAvgQueue << sampler.mNumQueues;
    return ar;
}

template <class T, class MatA, class MatB>
Archive& operator>>(Archive &ar, GibbsSampler<T, MatA, MatB> &sampler)
{
    ar >> sampler.mMatrix >> sampler.mAPMatrix >> sampler.mQueue >>
        sampler.mDomain >> sampler.mLambda >> sampler.mMaxGibbsMass >>
        sampler.mAnnealingTemp >> sampler.mNumRows >> sampler.mNumCols >>
        sampler.mBinSize >> sampler.mAvgQueue >> sampler.mNumQueues;
    return ar;
}

template <class T, class MatA, class MatB>
T* GibbsSampler<T, MatA, MatB>::impl()
{
    return static_cast<T*>(this);
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::processProposal(AtomicProposal *prop)
{
    unsigned r1 = 0, c1 = 0, r2 = 0, c2 = 0;
    switch (prop->type)
    {
        case 'B':
            birth(prop);
            break;
        case 'D':
            death(prop);
            break;
        case 'M':
            move(prop);
            break;
        case 'E':
            exchange(prop);
            break;
    }
}

// add an atom at pos, calculate mass either with an exponential distribution
// or with the gibbs mass calculation
template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::birth(AtomicProposal *prop)
{
    unsigned row = impl()->getRow(prop->birthPos);
    unsigned col = impl()->getCol(prop->birthPos);

    // calculate proposed mass
    float mass = 0.f;
    if (impl()->canUseGibbs(row, col))
    {
        AlphaParameters alpha = impl()->alphaParameters(row, col);
        mass = gibbsMass(alpha, &(prop->rng)).first;
    }
    else
    {
        mass = prop->rng.exponential(mLambda);
    }

    // accept mass as long as it's non-zero
    if (mass >= gaps::epsilon)
    {
        mDomain.cacheInsert(prop->birthPos, mass);
        mMatrix(row, col) += mass;
        impl()->updateAPMatrix(row, col, mass);
        mQueue.acceptBirth();
    }
    else
    {
        mQueue.rejectBirth();
    }
}

// automatically accept death, attempt a rebirth at the same position, using
// the original mass or the gibbs mass calculation
template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::death(AtomicProposal *prop)
{
    // calculate bin for this atom
    unsigned row = impl()->getRow(prop->atom1->pos);
    unsigned col = impl()->getCol(prop->atom1->pos);

    // kill off atom
    mMatrix(row, col) = gaps::max(mMatrix(row, col) - prop->atom1->mass, 0.f);
    impl()->updateAPMatrix(row, col, -1.f * prop->atom1->mass);

    // calculate rebirth mass
    float rebirthMass = prop->atom1->mass;
    AlphaParameters alpha = impl()->alphaParameters(row, col);
    if (impl()->canUseGibbs(row, col))
    {
        std::pair<float, bool> gMass = gibbsMass(alpha, &(prop->rng));
        if (gMass.second)
        {
            rebirthMass = gMass.first;
        }
    }

    // accept/reject rebirth
    float deltaLL = rebirthMass * (alpha.su - alpha.s * rebirthMass / 2.f);
    if (deltaLL * mAnnealingTemp >= std::log(prop->rng.uniform()))
    {
        prop->atom1->mass = rebirthMass;
        mMatrix(row, col) += rebirthMass;
        impl()->updateAPMatrix(row, col, rebirthMass);
        mQueue.rejectDeath();
    }
    else
    {
        mDomain.cacheErase(prop->atom1->pos);
        mQueue.acceptDeath();
    }
}

// move mass from src to dest in the atomic domain
template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::move(AtomicProposal *prop)
{
    unsigned r1 = impl()->getRow(prop->atom1->pos);
    unsigned c1 = impl()->getCol(prop->atom1->pos);
    unsigned r2 = impl()->getRow(prop->moveDest);
    unsigned c2 = impl()->getCol(prop->moveDest);

    GAPS_ASSERT(r1 != r2 || c1 != c2); // same bin handled during proposal

    float deltaLL = impl()->computeDeltaLL(r1, c1, -1.f * prop->atom1->mass,
        r2, c2, prop->atom1->mass);
    if (deltaLL * mAnnealingTemp > std::log(prop->rng.uniform()))
    {
        //mDomain.cacheMove(prop->atom1->pos, prop->moveDest, prop->atom1->mass);
        prop->atom1->pos = prop->moveDest; // temporary for vector domain

        mMatrix(r1, c1) = gaps::max(mMatrix(r1, c1) - prop->atom1->mass, 0.f);
        mMatrix(r2, c2) += prop->atom1->mass;

        impl()->updateAPMatrix(r1, c1, -1.f * prop->atom1->mass);
        impl()->updateAPMatrix(r2, c2, prop->atom1->mass);
    }
}

// exchange some amount of mass between two positions, note it is possible
// for one of the atoms to be deleted if it's mass becomes too small after
// the exchange
template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::exchange(AtomicProposal *prop)
{
    unsigned r1 = impl()->getRow(prop->atom1->pos);
    unsigned c1 = impl()->getCol(prop->atom1->pos);
    unsigned r2 = impl()->getRow(prop->atom2->pos);
    unsigned c2 = impl()->getCol(prop->atom2->pos);

    GAPS_ASSERT(r1 != r2 || c1 != c2); // same bin handled during proposal

    float m1 = prop->atom1->mass;
    float m2 = prop->atom2->mass;

    if (impl()->canUseGibbs(r1, c1, r2, c2))
    {
        AlphaParameters alpha = impl()->alphaParameters(r1, c1, r2, c2);
        std::pair<float, bool> gMass = gibbsMass(alpha, m1, m2, &(prop->rng));
        if (gMass.second)
        {
            acceptExchange(prop, gMass.first, r1, c1, r2, c2);
            return;
        }
    }

    float pUpper = gaps::p_gamma(m1 + m2, 2.f, 1.f / mLambda);
    float newMass = prop->rng.inverseGammaSample(0.f, pUpper, 2.f, 1.f / mLambda);

    float delta = m1 > m2 ? newMass - m1 : m2 - newMass; // change larger mass
    float pOldMass = 2.f * newMass > m1 + m2 ? gaps::max(m1, m2) : gaps::min(m1, m2);

    float pNew = gaps::d_gamma(newMass, 2.f, 1.f / mLambda);
    float pOld = gaps::d_gamma(pOldMass, 2.f, 1.f / mLambda);

    if (pOld == 0.f && pNew != 0.f) // special case
    {
        acceptExchange(prop, delta, r1, c1, r2, c2);
        return;
    }

    float deltaLL = impl()->computeDeltaLL(r1, c1, delta, r2, c2, -delta);
    float priorLL = (pOld == 0.f) ? 1.f : pOld / pNew;
    float u = std::log(prop->rng.uniform() * priorLL);
    if (u < deltaLL * mAnnealingTemp)
    {
        acceptExchange(prop, delta, r1, c1, r2, c2);
        return;
    }

    mQueue.rejectDeath();
}

// helper function for acceptExchange, used to conditionally update the mass
// at a single position, deleting the atom if the new mass is too small
template <class T, class MatA, class MatB>
bool GibbsSampler<T, MatA, MatB>::updateAtomMass(Atom *atom, float delta)
{
    if (atom->mass + delta < gaps::epsilon)
    {
        mDomain.cacheErase(atom->pos);
        mQueue.acceptDeath();
        return false;
    }
    atom->mass += delta;
    return true;
}

// helper function for exchange step, updates the atomic domain, matrix, and
// A*P matrix
template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::acceptExchange(AtomicProposal *prop,
float d1, unsigned r1, unsigned c1, unsigned r2, unsigned c2)
{
    // can make this assertion more formal if we customize p/q norm functions
    // to handle tail situations a certain way
    GAPS_ASSERT(prop->atom1->mass + d1 > gaps::epsilon);
    GAPS_ASSERT(prop->atom2->mass - d1 > gaps::epsilon);
    //d1 = gaps::max(-1.f * (prop->atom1->mass - epsilon), d1);

    float d2 = -1.f * d1;
    bool b1 = updateAtomMass(prop->atom1, d1);
    bool b2 = updateAtomMass(prop->atom2, d2);
    GAPS_ASSERT(b1 || b2);
    
    // delete entire atom if resize would make it too small
    if (!b1) { d1 = -1.f * prop->atom1->mass; }
    if (!b2) { d2 = -1.f * prop->atom2->mass; }

    if (b1 && b2)
    {
        mQueue.rejectDeath();
    }

    mMatrix(r1, c1) += d1;
    mMatrix(r2, c2) += d2;
    impl()->updateAPMatrix(r1, c1, d1);
    impl()->updateAPMatrix(r2, c2, d2);
}

template <class T, class MatA, class MatB>
std::pair<float, bool> GibbsSampler<T, MatA, MatB>::gibbsMass(AlphaParameters alpha,
GapsRng *rng)
{        
    alpha.s *= mAnnealingTemp;
    alpha.su *= mAnnealingTemp;

    if (alpha.s > gaps::epsilon)
    {
        float mean = (alpha.su - mLambda) / alpha.s;
        float sd = 1.f / std::sqrt(alpha.s);
        float pLower = gaps::p_norm(0.f, mean, sd);

        if (pLower < 1.f)
        {
            float m = rng->inverseNormSample(pLower, 1.f, mean, sd);
            float gMass = gaps::min(m, mMaxGibbsMass / mLambda);
            return std::pair<float, bool>(gMass, gMass >= gaps::epsilon);
        }
    }
    return std::pair<float, bool>(0.f, false);
}

template <class T, class MatA, class MatB>
std::pair<float, bool> GibbsSampler<T, MatA, MatB>::gibbsMass(AlphaParameters alpha,
float m1, float m2, GapsRng *rng)
{
    alpha.s *= mAnnealingTemp;
    alpha.su *= mAnnealingTemp;

    if (alpha.s > gaps::epsilon)
    {
        float mean = alpha.su / alpha.s; // lambda cancels out
        float sd = 1.f / std::sqrt(alpha.s);
        float pLower = gaps::p_norm(-m1, mean, sd);
        float pUpper = gaps::p_norm(m2, mean, sd);

        if (!(pLower >  0.95f || pUpper < 0.05f))
        {
            float delta = rng->inverseNormSample(pLower, pUpper, mean, sd);
            float gibbsMass = gaps::min(gaps::max(-m1, delta), m2); // conserve mass
            return std::pair<float, bool>(gibbsMass, true);
        }
    }
    return std::pair<float, bool>(0.f, false);
}

#endif
