extern "C" {
  #include "EXTERN.h"
  #include "perl.h"
  #include "XSUB.h"
  #undef seed
  #undef do_open
  #undef do_close
}

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include <interval_tree.h>

using namespace std;

#define do_open   Perl_do_open
#define do_close  Perl_do_close

class SV_ptr {
  SV *sv;
  public:
    SV_ptr() : sv(0) {}
    SV_ptr(SV *sv) : sv(sv) {
      if (sv) SvREFCNT_inc(sv);
    }
    SV_ptr(const SV_ptr &ptr) : sv(ptr.get()) {
      if (sv) SvREFCNT_inc(sv);
    }
    SV_ptr(const char* s){
    	sv = newSVpv(s,0);
    	if (sv) SvREFCNT_inc(sv);
    }
    virtual ~SV_ptr() {
      if (sv) SvREFCNT_dec(sv);
    }
    SV_ptr& operator=(SV_ptr ptr) {
      if (sv) SvREFCNT_dec(sv);
      sv = ptr.get();
      if (sv) SvREFCNT_inc(sv);
      return *this;
    }
    SV * get() {
      return sv; 
    }
    SV * get() const {
      return sv; 
    }
};

std::ostream& operator<<(std::ostream &out, SV_ptr value) {
  out << SvPV_nolen(value.get());
  return out;
}

class RemoveFunctor {
  SV *callback;
  public:
    RemoveFunctor(SV *callback_) : callback(callback_) {}
    bool operator()(SV_ptr value, long low, long high) const {
      // pass args into callback
      dSP;
      ENTER;
      SAVETMPS;
      PUSHMARK(SP);
      XPUSHs(value.get());
      XPUSHs(sv_2mortal(newSViv(low)));
      XPUSHs(sv_2mortal(newSViv(high))); // high + 1 deleted here #OD
      PUTBACK;

      // get result from callback and return
      I32 count = call_sv(callback, G_SCALAR);

      SPAGAIN;

      if (count < 1) {
        PUTBACK;
        FREETMPS;
        LEAVE;
        return false;
      }

      SV *retval_sv = POPs;
      bool retval = SvTRUE(retval_sv);

      PUTBACK;
      FREETMPS;
      LEAVE;
      return retval;
    }
};

typedef IntervalTree<SV_ptr,long> PerlIntervalTree;
typedef IntervalTree<SV_ptr,long>::Node PerlIntervalTree_Node;

MODULE = Set::IntervalTree PACKAGE = Set::IntervalTree

PerlIntervalTree *
PerlIntervalTree::new()

SV *
PerlIntervalTree::str()
  CODE:
    std::string str = THIS->str();
    const char *tree = str.c_str();
    RETVAL = newSVpv(tree, 0);
  OUTPUT:
    RETVAL

void
PerlIntervalTree::SaveTree(string filename)
    
void
PerlIntervalTree::LoadTree(string filename)

void
PerlIntervalTree::insert(SV *value, long low, long high)
  PROTOTYPE: $;$;$
  CODE: 
    if (high < low) Perl_croak(aTHX_ "Intervals must have positive width");
    SV_ptr ptr(value);
    THIS->insert(ptr, low, high);

AV *
PerlIntervalTree::remove(long low, long high, ...)
  CODE:
    if (high < low) Perl_croak(aTHX_ "Intervals must have positive width");
    RETVAL = newAV();
    sv_2mortal((SV*)RETVAL);

    if (items > 3) {
      SV *callback = ST(3); 
      RemoveFunctor remove_functor(callback);
      std::vector<SV_ptr> removed;
      THIS->remove(low, high, remove_functor, removed);

      for (std::vector<SV_ptr>::iterator
          i=removed.begin(); i!=removed.end(); ++i) 
      {
        SV *value = i->get();
        SvREFCNT_inc(value);
        av_push(RETVAL, value);
      }
    }
    else {
      std::vector<SV_ptr> removed; 
      THIS->remove(low, high, removed);

      for (std::vector<SV_ptr>::iterator
          i=removed.begin(); i!=removed.end(); ++i) 
      {
        SV *value = i->get();
        SvREFCNT_inc(value);
        av_push(RETVAL, value);
      }
    }
  OUTPUT:
    RETVAL

AV *
PerlIntervalTree::remove_window(long low, long high, ...)
  CODE:
    if (high < low) Perl_croak(aTHX_ "Intervals must have positive width");
    RETVAL = newAV();
    sv_2mortal((SV*)RETVAL);

    if (items > 3) {
      SV *callback = ST(3); 
      RemoveFunctor remove_functor(callback);
      std::vector<SV_ptr> removed;
      THIS->remove_window(low, high, remove_functor, removed);

      for (std::vector<SV_ptr>::iterator
          i=removed.begin(); i!=removed.end(); ++i) 
      {
        SV *value = i->get();
        SvREFCNT_inc(value);
        av_push(RETVAL, value);
      }
    }
    else {
      std::vector<SV_ptr> removed; 
      THIS->remove_window(low, high, removed);

      for (std::vector<SV_ptr>::iterator
          i=removed.begin(); i!=removed.end(); ++i) 
      {
        SV *value = i->get();
        SvREFCNT_inc(value);
        av_push(RETVAL, value);
      }
    }
  OUTPUT:
    RETVAL

AV *
PerlIntervalTree::fetch(long low, long high)
  PROTOTYPE: $;$
  CODE:
    if (high < low) Perl_croak(aTHX_ "Intervals must have positive width");
    RETVAL = newAV();
    sv_2mortal((SV*)RETVAL);
    std::vector<SV_ptr> intervals;
    THIS->fetch(low, high, intervals);
    for (size_t i=0; i<intervals.size(); i++) {
      SV *value = intervals[i].get();
      SvREFCNT_inc(value);
      av_push(RETVAL, value);
    }
  OUTPUT:
    RETVAL

AV *
PerlIntervalTree::fetch_window(long low, long high)
  PROTOTYPE: $;$
  CODE:
    if (high < low) Perl_croak(aTHX_ "Intervals must have positive width");
    RETVAL = newAV();
    sv_2mortal((SV*)RETVAL);
    std::vector<SV_ptr> intervals;
    THIS->fetch_window(low, high, intervals);
    for (size_t i=0; i<intervals.size(); i++) {
      SV *value = intervals[i].get();
      SvREFCNT_inc(value);
      av_push(RETVAL, value);
    }
  OUTPUT:
    RETVAL

void 
PerlIntervalTree::DESTROY()

MODULE = Set::IntervalTree PACKAGE = Set::IntervalTree::Node

PerlIntervalTree_Node *
PerlIntervalTree_Node::new()

int
PerlIntervalTree_Node::low()
  CODE:
    RETVAL = THIS->low();
  OUTPUT:
    RETVAL

int
PerlIntervalTree_Node::high()
  CODE:
    RETVAL = THIS->high();
  OUTPUT:
    RETVAL

SV *
PerlIntervalTree_Node::value()
  CODE:
    RETVAL = THIS->value().get();
  OUTPUT:
    RETVAL

void 
PerlIntervalTree_Node::DESTROY()

