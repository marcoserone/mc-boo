
\version "2.19.35"

\header {
  title = "Villagers' Chorus"
  subtitle = #(strftime "%d-%m-%Y" (localtime (current-time)))
  composer = "Borodin"
}

global = {
  \key a \major
  \time 4/4
}

#(set-global-staff-size 19)


soprano_high = \relative c'' {
  \global
fis4( cis'2) b8( a16 g)
fis4. fis8  b'4.( cis16 b
a4) fis e4.'( cis8
d e) cis( a) b4( fis,4)
cis'4. cis8 b ( cis) b(a)
e2(b'2) cis4.( e8 d  cis16 b16 a8 )b
fis2(cis')
b4( fis' e8 fis cis b)
a2 e'4.( d16 cis16
b2 fis4) cis'
e2. e4
e4(fis2.~fis1\fermata)


}
soprano_low = \relative c'' {
  \global
fis4,( cis'2) b8( a16 g)
fis4. fis8  b'4.( cis16 b
a4) fis e4.'( cis8
d e) cis( a) b4( fis,4)
cis'4. cis8 b ( cis) b(a)
e2(b'2) cis4.( e8 d  cis16 b16 a8 )b
fis2(cis')
b2.( e,4)
a2 a( 
b2 fis4) fis8( a)
b(cis) b( a) e4 e'
a,( cis8 e d cis16 b ) a8(b)
fis1\fermata


}

alto_high = \relative c' {
  \global

r4 r4 r4 r4 
fis8(g) a4 g fis
e2. e8 g
a8( g fis2.~fis2) fis2
fis2(g) fis  e(fis2.) a8(fis)
g(a) g(fis) e4 b'4
cis2( b8 a16 g) fis8(g) fis1\fermata
}
alto_low = \relative c' {
  \global

r4 r4 r4 r4 
fis8(g) a4 g fis
e2. e8 g
a8( g fis2.~fis2) fis4(e)
d2(e) fis4( d)  c( e
d b fis) fis
e e e e8(g) a( g fis2) fis8(g)
fis1\fermata
}

tenor = \relative c' {
  \global

  r4 cis8\pp cis16 cis~ cis8 cis16 cis~ cis8 cis16 cis~
  cis2 r
  r4 cis8 cis16 cis~ cis8 cis16 cis~ cis8 cis16 cis~
  cis2 r
  
  r4 cis8 cis cis cis cis cis
  b8 b b b16 b~\< b2
  
  r4\! cis8 cis16 cis~ cis8 cis16 cis~ cis8 cis16 cis~
  cis2 r
  r8. cis16 cis8 cis cis cis cis cis
  cis2 r
  
  r4 cis8 cis cis16 cis8. cis8 cis
  b b b b b2
  
  r1\fermata
  
  cis16 cis8. cis8 cis cis8. cis16 d8 e16 eis~
  eis4 r8 eis cis d cis cis16 d~
  d2 r4 d8 d16 d~
  d8 d d d16 cis~ cis8 cis cis d16 d~
  
  d2 r4 d8 d
  d8. d16 d8 d16 cis~ cis8 cis cis e16 d~
  d2 r\fermata
  
  
  e,16 e8. a8 b16 d~ d8 cis16 b~ b8 a16 gis~
  gis4 r8\fermata gis cis cis cis cis16 d~
  d2 r4 d16 d8.
  d8 d d d16 cis~ cis8 cis16 cis8. e16 d~
  
  d2 r4 d8 d
  d8 d d d cis16 cis cis8~ cis d16 d~
  d2 r
  
  d2. d8 d
  cis16 cis cis cis cis8 cis4. r4
  d2~ d8 d16 d~ d8 d16 cis~
  cis2 r
  
  b16 b8 b16~ b8 b
  cis r b8 a
  gis16 gis8 a16~ a8 gis
  a2
  
  r1
  
  cis8 cis cis cis cis cis cis cis16 cis~
  cis4 r8 cis cis d cis cis16 d~
  d2 r4 d8 d
  d8 d d d16 cis~ cis8 cis cis d16 d~
  
  d2 r4 d8 d
  d8 d d d cis16 cis cis8 r cis16 d~
  d2 r
 
  b8 b b8. b16 b8 b16 b~ b8 cis
  b16 cis8 cis16 r8 cis cis16 b8 a16 r8 gis16 a
  b8 b b a16 b~ b4 r8\fermata b16 b
  a8 a a a16 a~ a4 r8\fermata a16 a
  b8 b b cis16 b~ b4 r8\fermata b16 cis
  
  d2. d8 d
  cis16 cis cis cis cis8 cis4. r4
  d2~ d8 d16 d~ d8 d16 cis~
  cis2 r
  
  b16\mp\> b8 b16~ b8 b
  cis r b8 a
  gis16 gis8 a16~ a8 gis
  a2\!
}

bass = \relative c' {
  \global
  R1*6
  
  
  \override NoteHead.style = #'cross
  a4 \pp  r a r
  a r a r
  a r a r
  a r a r
  
  a r a r
  a a r2
  \revert NoteHead.style
  
  r1
  
  a16\mf a8. a8 a a8. a16 a8 a16 cis~
  cis4 r8 cis gis gis gis gis16 fis~
  fis2 r4 fis8 gis16 a~
  a8 a a gis16 a~ a8 a a gis16 a~
  
  a2 r4 a8 gis
  a8. a16 a8 gis16 a~ a8 a a a16 a~
  a2 r2
  
  r1
  r4. gis8 gis gis gis gis16 fis~
  fis2 r4 fis16 gis8.
  a8 a a gis16 a~ a8 a16 a8. gis16 a~
  
  a2 r4 a8 gis
  a8 a a8 gis a16 a a8~ a8 a16 a~
  a2 r2
  
  
  d2.\pp\< d8 d
  cis16 cis cis cis cis8 cis4.\! r4
  d2~\pp\< d8 d16 d~ d8 d16 cis~
  cis2 r\!
  
  g16\mp g8 g16~ g8 g
  fis r fis fis
  e16 e8 e16~ e8 e
  d2
  
  r1
  
  a'8 a a a a a a a16 gis~
  gis4 r8 gis gis gis gis gis16 fis~
  fis2 r4 fis8 gis
  a8 a a gis16 a~ a8 a a gis16 a~
  
  a2 r4 a8 gis
  a8 a a8 gis a16 a a8 r a16 a~
  a2 r2
  
  
  e8\f e e8. e16 eis8 eis16 eis~ eis8 eis
  fis16 fis8 fis16 r8 fis e16 e8 e16 r8 e16 e
  dis8 dis dis dis16 dis~ dis4 r8 dis16\> dis
  e8 e e e16 e~ e4 r8 e16 fis
  gis8 gis gis gis16 gis~ gis4 r8 b16 cis
  
  d2.\pp\< d8 d
  cis16 cis cis cis cis8 cis4.\! r4
  d2~\pp\< d8 d16 d~ d8 d16 cis~
  cis2 r\!
  
  g16 g8 g16~ g8 g
  fis r fis fis
  e16 e8 e16~ e8 e
  d2
}

verseOne = \lyricmode {
  Ground Con -- trol to Ma -- jor Tom
  Ground Con -- trol to Ma -- jor Tom
  take your pro -- tein pills and put your hel -- met on
  Ground Con -- trol to Ma -- jor Tom
  com -- men -- cing count -- down, en -- gines on
  check ig -- ni -- tion and may God's love be with you
}

verseTwo = \lyricmode {
  this is Ground Con -- trol to Ma -- jor Tom
  you've real -- ly made the grade
  and the pa -- pers want to know whose shirts you wear
  now it's time to leave the cap -- sule if you dare
}

verseThree = \lyricmode {
  I'm step -- ping through the door
  and I'm floa -- ting in a most pe -- cu -- liar way
  and the stars look ve -- ry dif -- fe -- rent to -- day
}

bridge = \lyricmode {
  here am I sit -- ting in a tin can
  far a -- bove the world
  pla -- net Earth is blue
  and there's no -- thing I can do
}

verseFour = \lyricmode {
  though I'm past one hun -- dred thou -- sand miles
  I'm fee -- ling ve -- ry still
  and I think my space -- ship knows which way to go
  tell my wife I love her ve -- ry much she knows
}

verseFive = \lyricmode {
  Ground Con -- trol to Ma -- jor Tom
  your cir -- cuit's dead, there's some -- thing wrong,
  can you hear me, Ma -- jor Tom?
  can you hear me, Ma -- jor Tom?
  can you hear me, Ma -- jor Tom?
  can you
}

sopranoVerse = \lyricmode {
  \verseOne
  \verseTwo
  \verseThree
  \bridge
  \verseFour
  \verseFive
  \bridge
}

altoVerse = \lyricmode {
  \verseOne
  \verseTwo
  \verseThree
  \bridge
  \verseFour
  \verseFive
  \bridge
}

tenorVerse = \lyricmode {
  \verseOne
  \verseTwo
  This is Ma -- jor Tom to Ground Con -- trol
  \verseThree
  \bridge
  \verseFour
  \verseFive
  \bridge
}

bassVerse = \lyricmode {
  ten nine eight seven six five four three two one lift off
  \verseTwo
  \verseThree
  \bridge
  \verseFour
  \verseFive
  \bridge
}

chordsPart = \new ChordNames \chordNames

choirPart = \new ChoirStaff <<
  \new Staff = "sa" \with {
    instrumentName = \markup \center-column { "Sopran" "Alt" }
  } <<
    \new Voice = "soprano" { \voiceOne \soprano }
    \new Voice = "alto" { \voiceTwo \alto }
  >>
  \new Lyrics \with {
    alignAboveContext = "sa"
    \override VerticalAxisGroup #'staff-affinity = #DOWN
  } \lyricsto "soprano" \sopranoVerse
  \new Lyrics \lyricsto "alto" \altoVerse
  \new Staff = "tb" \with {
    instrumentName = \markup \center-column { "Tenor" "Bass" }
  } <<
    \clef bass
    \new Voice = "tenor" { \voiceOne \tenor }
    \new Voice = "bass" { \voiceTwo \bass }
  >>
  \new Lyrics \with {
    alignAboveContext = "tb"
    \override VerticalAxisGroup #'staff-affinity = #DOWN
  } \lyricsto "tenor" \tenorVerse
  \new Lyrics \lyricsto "bass" \bassVerse
>>

\score {
  <<
    \chordsPart
    \choirPart
  >>
  \layout { }
  \midi {
    \tempo 4=60
  }
}
