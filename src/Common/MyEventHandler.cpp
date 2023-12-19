// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
#pragma warning(disable : 4786)
#pragma warning(disable : 4503)
#endif

#include "MyEventHandler.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
MyEventHandler::MyEventHandler()
  : ClpEventHandler()
{
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
MyEventHandler::MyEventHandler(const MyEventHandler &rhs)
  : ClpEventHandler(rhs)
{
}

// Constructor with pointer to model
MyEventHandler::MyEventHandler(ClpSimplex *model)
  : ClpEventHandler(model)
{
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
MyEventHandler::~MyEventHandler()
{
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
MyEventHandler &
MyEventHandler::operator=(const MyEventHandler &rhs)
{
  if (this != &rhs) {
    ClpEventHandler::operator=(rhs);
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpEventHandler *MyEventHandler::clone() const
{
  return new MyEventHandler(*this);
}

int MyEventHandler::event(Event whichEvent)
{
    // This is how one can get some progress information at the end of each iteration.
    if (whichEvent == endOfIteration) {
        int numIter = model_->getIterationCount();
        int varIn = model_->sequenceIn();
        int varOut = model_->sequenceOut();
        int pivotRow = model_->pivotRow();
        int dirIn = model_->directionIn();
        int dirOut = model_->directionIn();
        int numStruct = model_->getNumCols();
        int numArtificial = model_->getNumRows();
        int status = model_->isProvenPrimalInfeasible();
        double *rowStruct = new double[numStruct];
        double *rowArtificial = new double[numArtificial]; 
        model_->getBInvARow(pivotRow, rowStruct, rowArtificial);
        std::cout << 
        "++++++++++++++++++++++++++++" <<
        "Iter: " << numIter << "\n" <<
        "Variable In: " << varIn << ", direction " << dirIn << "\n" <<
        "Variable Out: " << varOut << ", direction " << dirOut << "\n" <<
        "Pivot Row: " << pivotRow << "\n";

    }

    return -1;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
