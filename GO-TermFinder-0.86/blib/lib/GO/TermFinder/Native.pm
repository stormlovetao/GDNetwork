# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.35
#
# Don't modify this file, modify the SWIG interface instead.

package GO::TermFinder::Native;
require Exporter;
require DynaLoader;
@ISA = qw(Exporter DynaLoader);
package GO::TermFinder::Nativec;
bootstrap GO::TermFinder::Native;
package GO::TermFinder::Native;
@EXPORT = qw( );

# ---------- BASE METHODS -------------

package GO::TermFinder::Native;

sub TIEHASH {
    my ($classname,$obj) = @_;
    return bless $obj, $classname;
}

sub CLEAR { }

sub FIRSTKEY { }

sub NEXTKEY { }

sub FETCH {
    my ($self,$field) = @_;
    my $member_func = "swig_${field}_get";
    $self->$member_func();
}

sub STORE {
    my ($self,$field,$newval) = @_;
    my $member_func = "swig_${field}_set";
    $self->$member_func($newval);
}

sub this {
    my $ptr = shift;
    return tied(%$ptr);
}


# ------- FUNCTION WRAPPERS --------

package GO::TermFinder::Native;


############# Class : GO::TermFinder::Native::Distributions ##############

package GO::TermFinder::Native::Distributions;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( GO::TermFinder::Native );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = GO::TermFinder::Nativec::new_Distributions(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        GO::TermFinder::Nativec::delete_Distributions($self);
        delete $OWNER{$self};
    }
}

*pValueByHypergeometric = *GO::TermFinder::Nativec::Distributions_pValueByHypergeometric;
*hypergeometric = *GO::TermFinder::Nativec::Distributions_hypergeometric;
*logNCr = *GO::TermFinder::Nativec::Distributions_logNCr;
*logFactorial = *GO::TermFinder::Nativec::Distributions_logFactorial;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


# ------- VARIABLE STUBS --------

package GO::TermFinder::Native;

1;
